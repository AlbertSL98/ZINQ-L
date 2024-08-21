#' @import lme4
#' @import dplyr
#' @import quantreg
#' @import Matrix
#' @import MASS

score_func = function(tau, u){
  tau - ifelse(u < 0, 1, 0)
}


QRank_test = function(X1m,Z, study, y, ybin, taus = c(0.1, 0.25, 0.5, 0.75, 0.9)){
  qfit = QRank(gene = y, snp = X1m, cov = Z, tau = taus)
  pval = qfit$composite.pvalue
  return(pval)
}

logit_mixed_test = function(X1m, Z, study, y, ybin){
  logit_model = try(glmer(ybin ~ -1 + X1m + Z + (1|study),
                          family = binomial(link = "logit")), silent = T)

  if ('try-error' %in% class(logit_model)){
    logit_pval = NA
    names(logit_pval) = 'logit_Wald'
    logit_pval_LRT = NA
    names(logit_pval_LRT) = 'logit_LRT'
  } else{
    logit_summ = summary(logit_model)
    logit_pval = logit_summ$coefficients['X1m', 'Pr(>|z|)']
    names(logit_pval) = 'logit_Wald'

    logit_model_null = try(glmer(ybin ~ -1 + Z + (1|study),
                                 family = binomial(link = "logit")), silent = T)
    logit_pval_LRT = anova(logit_model, logit_model_null)$`Pr(>Chisq)`[2]
    names(logit_pval_LRT) = 'logit_LRT'
  }
  pval = c(logit_pval, logit_pval_LRT)
  return(pval)
}


minP_combine_test = function(pval,
                             taus = c(0.1, 0.25, 0.5, 0.75, 0.9),
                             Sigma.hat,
                             M = 10000, ignore_logit=0){
  l = length(pval)
  width = length(taus)
  df = 1
  if(ignore_logit==1){
	pval = pval[2:l]
  }
  t.obs = min(pval)
  qmin.quantile = rep(qchisq(1-t.obs, df= df), width)

  eSig = eigen(Sigma.hat,symmetric = T)
  eSig.value = eSig$values
  if (any(eSig.value < -1e-10)){
    eSig.vector = eSig$vectors[,which(eSig.value>1e-10)]
    Lambda = as.matrix(diag(eSig.value[which(eSig.value>1e-10)]))
    Sigma.hat = eSig.vector %*% Lambda %*% t(eSig.vector)
  }
  beta.sim = MASS::mvrnorm(n=M, mu=rep(0, width*df), Sigma=Sigma.hat)
  prob.quantile = mean(apply(beta.sim, 1, function(z){
    obs.sim = z^2/diag(Sigma.hat)
    all(obs.sim < qmin.quantile)
  }) )
  if(ignore_logit==0){
	min_pval = 1 - (1-t.obs)*prob.quantile
  } else{
    min_pval = 1 - prob.quantile
  }
  return(min_pval)

}

Truncated_Cauchy_Combination_zerorate_equal <-function(zerorate, pvalues, w = NULL, delta = 0.01, ignore_logit=0){
  l = length(pvalues)
  if(is.null(w)){
  	w = rep(-1,l)
  	w[1] = zerorate
  	w[2:l] = rep( (1-zerorate) / (l-1), l-1 )
  	w = w /sum(w)
  }

  if(ignore_logit==1){
  	w = w[2:l]
  	pvalues = pvalues[2:l]
  	w = w / sum(w)
  }

  is.small <- ifelse( (pvalues < 1e-15)&(!is.na(pvalues)), TRUE, FALSE )
  is.big <- ifelse( (pvalues > 1-delta)&(!is.na(pvalues)), TRUE, FALSE )
  pvalues <- as.matrix(pvalues)
  pvalues[is.big] <- 1-delta
  pvalues[!is.small] <- tan((0.5-pvalues[!is.small])*pi)
  pvalues[is.small] <- 1/pvalues[is.small]/pi
  cct.stat <- sum(w*(pvalues),na.rm = F)
  pcc <- pcauchy(cct.stat,lower.tail = F)
  return(pcc)
}



#' @export
ZINQL_fit = function(y, formula, formula.logistics=NA, meta, C, taus=c(0.1, 0.25, 0.5, 0.75, 0.9), method='Both', n.positive.cut=5, seed=2024){
  set.seed(2024)
  #
  #parse method
  #method should be one of 'Both', 'MinP','Cauchy'
  if(method=='Both'){
    run.minP = 1
    run.cauchy = 1
    n.p.values = 2
    re.col.n = c('ZINQL_MinP','ZINQL_Cauchy')
    model.type='Both'
  } else if(method=='MinP'){
    run.minP = 1
    run.cauchy = 0
    n.p.values = 1
    re.col.n = c('ZINQL_MinP')
    model.type='ZINQL_MinP'
  } else if(method=='Cauchy'){
    run.minP = 0
    run.cauchy = 1
    n.p.values = 1
    re.col.n = c('ZINQL_Cauchy')
    model.type='ZINQL_Cauchy'
  } else{
    print('Method should be one of Both, MinP or Cauchy')
    result.full[['Final_P_value']] = NA
    result.full[['model']] ='None'
    result.full[['Intermediate_P_value']] = NA
    return(result.full)
  }
  #
  result.full = list()
  #
  #formula = y ~ X + Z + (1|study)
  #get variable names for fixed effect and random effect
  fixed_effects <- attr(terms(formula), "term.labels")
  fixed_effects <- fixed_effects[!grepl("\\|", fixed_effects)]
  fixed_effects = c(fixed_effects)
  #
  # Extract random effect grouping factor
  random_effects <- attr(terms(formula), "term.labels")
  random_effects <- random_effects[grepl("\\|", random_effects)]
  random_effects <- gsub(".*\\|", "", random_effects)
  random = trimws(random_effects)
  #
  # get variables other than C to adjust for
  Z = setdiff(fixed_effects, C)
  X = C
  #
  #method should be one of 'Both', 'MinP','Cauchy'
  #ybin is the binary for y before dither
  ybin = (y > 0) * 1
  #

  yd = dither(y, type = "right", value = 1)

  X1m = as.matrix( meta[,X] )
  colnames(X1m) = X
  Z1m = as.matrix( meta[,Z] )
  colnames(Z1m) = Z
  study = meta[,random]
  #
  #decide whether to add intercept
  formula_str = Reduce(paste, deparse(formula))
  has_no_intercept = grepl("-1", formula_str)
  if(!(has_no_intercept)){
    Z1m = cbind(1, Z1m)
    colnames(Z1m) = c( 'intercept',Z )
  }
  #first select variables globally. If one variable, in addition to intercept, has all unique value. We should drop it

  global.Z.keep = c()
  for(v.tmp in colnames(Z1m)){
    if(v.tmp == 'intercept'){
      global.Z.keep = c(global.Z.keep, 'intercept')
    } else{
      z.unique = unique(Z1m[,v.tmp])
      if(length(z.unique)>1){
        global.Z.keep=c(global.Z.keep,v.tmp)
      } else{
        print(paste0('The covariate ',v.tmp, ' has unique value, thus it is excluded in the analysis!'))
        }
    }
  }

  if(length(global.Z.keep)==0){
    print('There is no covariate left in addition to the variable of interest, and the user does not allow the existence of intercept, so the package will return NA. The user may at least include an intercept term!')
    result.full[['Final_P_value']] = rep(NA,run.minP+run.cauchy)
    result.full[['model']] ='None'
    result.full[['Intermediate_P_value']] = rep(NA,length(taus)+1)
    return(result.full)
  }

  x.unique = unique(meta[,X])
  if(length(x.unique)==1){
    print( paste0('The variable of interest: ', X ,', has only unique value, thus the model will return NA.') )
    result.full[['Final_P_value']] = rep(NA,run.minP+run.cauchy)
    result.full[['model']] ='None'
    result.full[['Intermediate_P_value']] = rep(NA,length(taus)+1)
    return(result.full)
  }

  Z1m = as.matrix(Z1m[,global.Z.keep])
  colnames(Z1m) = global.Z.keep
  #
  m = length(unique(study))
  N = length(X1m)
  ni = as.vector(table(study))
  L = sum(ni * (ni - 1))
  #q = ncol(Z1m)
  #
  nonzero_idx = which(ybin > 0)
  #print(length(nonzero_idx))

  Intermediate.col.name = c('Logistics_LRT')
  Intermediate.col.name = c(Intermediate.col.name, paste0('Quantile_',taus))


	#here we count number of zeros
	n.positive = length(which(y>0))
	#these two are indicators of specific procedures
	rely.logistics = 0
	fill.in.one.quantile = 0

	if(n.positive < n.positive.cut){
  	print( paste0('There are only ',n.positive, ' positive values in y, while ZINQ-L requires at least ', n.positive.cut, ' positive values. The package will return NA.') )
  	result.full[['Final_P_value']] = rep(NA,run.minP+run.cauchy)
  	result.full[['model']] ='None'
  	result.full[['Intermediate_P_value']] = rep(NA,length(taus)+1)
  	return(result.full)
	}

	#this part is to select variables for quantile regression, since it is possible that Z has unique value on POSITIVE part of y
	if(rely.logistics==0){
  	#now select Z. This is because sometimes positive part of y's Z may be unique, like for positive part of y, all gender are Female, which will result in failure of rq
	  quantile.Z.keep = c()
	  non.zero.index = which(y>0)
	  for(v.tmp in colnames(Z1m)){
	    if(v.tmp == 'intercept'){
	      quantile.Z.keep = c(quantile.Z.keep, 'intercept')
	    } else{
	      feature.value.unique = unique(meta[non.zero.index,v.tmp])
	      if(length(feature.value.unique)>1){
	        quantile.Z.keep=c(quantile.Z.keep, v.tmp)
	      } else{
	        print(paste0('The covariate ',v.tmp, ' has unique value for the POSITIVE part of y, thus it is excluded in the Quantile mixed effect Regression!'))
	      }
	    }
	  }
  	#
  	if(length(quantile.Z.keep)==0){
    	print('For the Quantile mixed effect Regression part, there is no covariate left in addition to the variable of interest, and the user does not allow the existence of intercept, so the test will be based on Logistics mixed effect Regression only. The user should at least include the intercept term, or choose a different set of covariates!')
    	rely.logistics = 1
  	}

  	Z1m.quantile = as.matrix( Z1m[,quantile.Z.keep] )
    colnames(Z1m.quantile) = quantile.Z.keep
	}

	x.quantile.unique = unique(meta[non.zero.index,X])
	if(length(x.quantile.unique)==1){
	  fill.in.one.quantile = 1
	  print( paste0('Variable of interest: ',X, ', has unique value on POSITIVE part of y, thus all the P-values for Quantile mixed effect Regression will be one.') )
	  if(run.minP==1){
	    print('The ZINQ-L MinP will not work in this situation, because the correlation between test statistics at different quantile levels cannot be properly estimated!')
	  }
	}

	X1m_tilde = X1m * ybin
	Z_tilde = Z1m.quantile * ybin

	print('Conducting Logistics mixed effect Regression.')
  #we provide the option that logistics regression formula is different from the quantile formula
  if(is.na(formula.logistics)){
  	## GLMM - logit reg
  	logit_model = try(glmer(ybin ~ -1 + X1m + Z1m + (1|study),
  							family = binomial(link = "logit")), silent = T)
  } else{
    #get variable names for fixed effect and random effect
    fixed_effects.logis <- attr(terms(formula.logistics), "term.labels")
    fixed_effects.logis <- fixed_effects.logis[!grepl("\\|", fixed_effects.logis)]
    fixed_effects.logis = c(fixed_effects.logis)
    # Extract random effect grouping factor
    random_effects.logis <- attr(terms(formula.logistics), "term.labels")
    random_effects.logis <- random_effects.logis[grepl("\\|", random_effects.logis)]
    random_effects.logis <- gsub(".*\\|", "", random_effects.logis)
    random.logis = trimws(random_effects.logis)
    # get variables other than C to adjust for
    Z.logis = setdiff(fixed_effects.logis, C)
    X.logis = C
    #for ligistics only
    X1m.logis = as.matrix( meta[,X.logis] )
    colnames(X1m.logis) = X.logis
    Z1m.logis = as.matrix( meta[,Z.logis] )
    colnames(Z1m.logis) = Z.logis
    study.logis = meta[,random.logis]
	#
	#decide whether to add intercept
    formula_str.logis = deparse(formula.logistics)
    has_no_intercept.logis = grepl("-1", formula_str.logis)
    if(!(has_no_intercept.logis)){
      Z1m.logis = cbind(1, Z1m.logis)
      colnames(X1m.logis) = c( 'intercept',Z.logis )
    }
    #fit logis
    ## GLMM - logit reg
    logit_model = try(glmer(ybin ~ -1 + X1m.logis + Z1m.logis + (1|study.logis),
                            family = binomial(link = "logit")), silent = T)
  }

	if ('try-error' %in% class(logit_model)){
	  logit_pval = NA
	  names(logit_pval) = 'logit_Wald'
	  logit_pval_LRT = NA
	  names(logit_pval_LRT) = 'logit_LRT'
	} else{
	  logit_summ = summary(logit_model)
	  logit_pval = logit_summ$coefficients['X1m', 'Pr(>|z|)']
	  names(logit_pval) = 'logit_Wald'

	  logit_model_null = try(glmer(ybin ~ -1 + Z1m + (1|study),
								   family = binomial(link = "logit")), silent = T)
	  logit_pval_LRT = anova(logit_model, logit_model_null)$`Pr(>Chisq)`[2]
	  names(logit_pval_LRT) = 'logit_LRT'
	}

	if(rely.logistics==1){
	  tmp = rep(logit_pval_LRT,n.p.values)
	  names(tmp) = re.col.n
	  result.full[['Final_P_value']] = tmp
	  result.full[['model']] ='Logistics'
	  tmp = c(logit_pval_LRT,rep(NA,length(Intermediate.col.name)-1))
	  names(tmp) = Intermediate.col.name
	  result.full[['Intermediate_P_value']] = tmp
	  return(result.full)
	}

	#whether ignore logit part
	ignore_logit = 0
	if( min(ybin)>0 ) {
		print('All the ys are positive, thus we manually set P-value for Logistics mixed effect Regression as 1!')
	  logit_pval_LRT = 1
	  names(logit_pval_LRT) = 'logit_LRT'
	}else if(is.na(logit_pval_LRT)){
	  print('Logistics mixed effect Regression returns NA P-value, thus we manually set the P-value as 1!')
	  logit_pval_LRT = 1
	  names(logit_pval_LRT) = 'logit_LRT'
	}

	print('Conducting Quantile mixed effect Regression.')
  X1m_star = try(X1m_tilde - Z_tilde %*% solve(crossprod(Z_tilde), crossprod(Z_tilde, X1m_tilde)), silent = T)
  if ('try-error' %in% class(X1m_star)){
    ZtZ= crossprod(Z_tilde)
    eZ = eigen(ZtZ,symmetric = T)
    eZ.value = eZ$values
    if (any(abs(eZ.value) < 1e-10)){
      n = length(which(abs(eZ.value) > 1e-10))
      eZ.vector = eZ$vectors[,which(abs(eZ.value)>1e-10), drop = F]
      if (n == 1){
        Lambda_inv = matrix(eZ.value[which(abs(eZ.value)>1e-10)]^(-1), 1, 1)
      } else{
        Lambda_inv = as.matrix(diag(eZ.value[which(abs(eZ.value)>1e-10)]^(-1)))
      }
      ZtZ_inv = eZ.vector %*% Lambda_inv %*% t(eZ.vector)
    }
    X1m_star = X1m_tilde - Z_tilde %*% ZtZ_inv %*% crossprod(Z_tilde, X1m_tilde)
  }
  #
  #df = data.frame(yd, ybin, X1m = X1m, Z1m.quantile, study, idx = 1:length(yd)) %>% filter(ybin > 0)
  df = data.frame(yd, ybin, X1m = X1m, Z1m.quantile, study, idx = 1:length(yd))
  df = df[which(df$ybin > 0),]

  nonzero_ni = as.vector(table(df$study))
  nonzero_L = sum(nonzero_ni * (nonzero_ni-1))
  #
  qformula = as.formula(paste0("yd ~ -1 + ", paste0(colnames(Z1m.quantile), collapse = '+')))
  qmf = model.frame(qformula, data= df)

  #if the quantile regression cannot be fit successfully. We also return P-values from logistics regression
  #print('Fitting the quantile regression.')
  qfit = try(rq(qformula, tau = taus, data = df), silent = T)
  if ('try-error' %in% class(qfit)){
    print('Quantile mixed effect Regression cannot be successfully fitted on POSITIVE part of y, thus all the P-values for Quantile mixed effect Regression will be manually set as 1!')
    fill.in.one.quantile = 1
    if(run.minP==1){
      print('The ZINQ-L MinP will not work in this situation, because the correlation between test statistics at different quantile levels cannot be properly estimated!')
    }
  }

  if(fill.in.one.quantile == 0){
    #qfit = rq(qformula, tau = taus, data = df)
    gamma = coef(qfit)
    #print(gamma)
    u = df$yd - predict(qfit)
    S = sapply(1:length(taus), function(i){
      1/sqrt(N) * sum(X1m_star[nonzero_idx] * score_func(taus[i],u[,i]))
    })
    #
    ## combine multiple quantiles
    width = length(taus)
    V0 = matrix(0, ncol=width, nrow=width)
    for (kk in 1:(width-1)){
      for (ll in (kk+1):width){
        V0[kk, ll] = min(taus[kk], taus[ll]) - taus[kk]*taus[ll]
      }
    }
    V0 = V0 + t(V0) + diag(taus*(1 - taus))
    V0 = V0 * sum(X1m_star[nonzero_idx]^2) / N

    Q_dep_result = lapply(1:length(taus), function(t){
      tau = taus[t]
      diag_ele = tau*(1-tau)* sum(X1m_star[nonzero_idx]^2)
      #
      off_diag_ele_result = lapply(unique(df$study), function(i){
        idx = df$idx[which(df$study == i)]
        #
        Xlist = lapply(idx, function(id) X1m_star[id,])
        Xsum = Reduce('+', Xlist)
        t1 = Xsum %*% t(Xsum)  - t(X1m_star[idx,]) %*% X1m_star[idx,]
        #
        Xlist_neg_res = lapply(idx, function(id) X1m_star[id,] * (u[which(df$idx == id),t] <0))
        X1m_neg_res = do.call(rbind, Xlist_neg_res)
        Xsum_neg_res = Reduce('+', Xlist_neg_res)
        t2 = Xsum_neg_res %*% t(Xsum)
        t2m = lapply(1:length(Xlist), function(l){
          Xlist_neg_res[[l]] %*% t(Xlist[[l]])
        }) %>% Reduce('+', .)
        t2 = t2 - t2m

        t3 = Xsum_neg_res %*% t(Xsum_neg_res)  - t(X1m_neg_res) %*% X1m_neg_res

        ans = list(t1 = t1,
                   t2 = t2,
                   t = tau^2 * t1 - 2 * tau* t2 + t3,
                   Xlist_neg_res = Xlist_neg_res
                   )
        return(ans)
      })

      off_diag_ele = lapply(off_diag_ele_result, function(i) i$t) %>% Reduce('+', .)
      t1 = lapply(off_diag_ele_result, function(i) i$t1) %>% Reduce('+', .)
      t2 = lapply(off_diag_ele_result, function(i) i$t2) %>% Reduce('+', .)

      Xlist_neg_res = lapply(off_diag_ele_result, function(i) i$Xlist_neg_res)

      Q_dep = as.numeric((diag_ele + off_diag_ele)/N)

      ans = list(t1 = t1,
                 t2 = t2,
                 Q_dep = Q_dep,
                 Xlist_neg_res = Xlist_neg_res)
    })

    Q_dep = sapply(Q_dep_result, function(q) q$Q_dep)
    Xlist_neg_res = lapply(Q_dep_result, function(q) q$Xlist_neg_res)
    t1 = sapply(Q_dep_result, function(q) q$t1)
    t2 = sapply(Q_dep_result, function(q) q$t2)

    ## combine multiple quantiles
    width = length(taus)
    V = matrix(0, width, width)
    for(i in 1:(width-1)){
      for(j in (i+1):width){
        tau_i = taus[i]
        tau_j = taus[j]
        d0 = (min(tau_i, tau_j) - tau_i * tau_j) * sum(X1m_star[nonzero_idx]^2)
        d1 = tau_i * tau_j * t1[1]
        d2 = -tau_j * t2[i]-tau_i * t2[j]
        li = Xlist_neg_res[[i]]
        lj = Xlist_neg_res[[j]]
        d3 = lapply(1:length(li), function(k){
          tmp_i = Reduce('+', li[[k]])
          tmp_j = Reduce('+', lj[[k]])
          d3 = tmp_i %*% t(tmp_j)  - t(do.call(rbind, li[[k]])) %*% do.call(rbind, lj[[k]])
        }) %>% Reduce('+', .)
        V[i, j] = (d0 + d1 + d2 + d3)/N
      }
    }
    V = V + t(V) + diag(Q_dep)
    #
    Ttau_dep = S^2/Q_dep
    pval_dep = sapply(Ttau_dep, function(tt){
      pchisq(tt, df = 1, lower.tail = F)
    })
    #names(pval_dep)= paste0("dep:", taus)
  }

  if(fill.in.one.quantile==1){
    pval_dep = rep(1,length(taus))
    V = diag(0,length(taus))
  }

  pval_dep = ifelse(is.nan(pval_dep),1,pval_dep)

  if(run.minP == 1){
	  pval_dep_minP = minP_combine_test(pval = c(logit_pval_LRT, pval_dep),
										taus = taus,
										Sigma.hat = V,
										M = 10000, ignore_logit=ignore_logit)
	  names(pval_dep_minP) = 'ZINQL_MinP'
  }

  zerorate = ifelse(is.na(logit_pval_LRT), 0, mean(ybin == 0))
  if(run.cauchy == 1){
    pval_dep_cauchy_trunc_zero_equal = Truncated_Cauchy_Combination_zerorate_equal(zerorate = zerorate, pvalues = c(logit_pval_LRT, pval_dep), ignore_logit=ignore_logit)
	names(pval_dep_cauchy_trunc_zero_equal) = "ZINQL_Cauchy"
  }

  if(method=='Both'){
  	res = c(pval_dep_minP, pval_dep_cauchy_trunc_zero_equal)
  	result.full[['Final_P_value']] = res
  	result.full[['model']] = model.type
  	tmp = c(logit_pval_LRT,pval_dep)
  	names(tmp) = Intermediate.col.name
  	result.full[['Intermediate_P_value']] = tmp
  	return(result.full)
	} else if(method=='MinP'){
	  res = c(pval_dep_minP)
	  result.full[['Final_P_value']] = res
	  result.full[['model']] = model.type
	  tmp = c(logit_pval_LRT,pval_dep)
	  names(tmp) = Intermediate.col.name
	  result.full[['Intermediate_P_value']] = tmp
	  return(result.full)
	} else if(method=='Cauchy'){
	  res = c(pval_dep_cauchy_trunc_zero_equal)
	  result.full[['Final_P_value']] = res
	  result.full[['model']] = model.type
	  tmp = c(logit_pval_LRT,pval_dep)
	  names(tmp) = Intermediate.col.name
	  result.full[['Intermediate_P_value']] = tmp
	  return(result.full)
	} else{ print('Method should be one of Both, MinP or Cauchy')
	}
}
