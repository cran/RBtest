#' Test of missing data mechanism using complete data
#'
#' This function tests the missing completely at random (MCAR) vs missing at random (MAR) by using the complete variables only.
#' @usage RBtest(data)
#' @param data Dataset with at least one complete variable. The variables could be either continuous, categorical or a mix of both.
#'
#' @return A list of the following elements:
#'  \itemize{
#'   \item \code{abs.nbrMD} The absolute number of missing data per variable.
#'   \item \code{rel.nbrMD} The percentage of missing data per variable.
#'   \item \code{type} Vector of the same length than the number of variables of the dataset, where '0' is for variables with MCAR data, '1' is for variables with MAR data and '-1' is for complete variables.
#'   }
#'
#' @examples
#'
#' set.seed(60)
#' n<-100 # sample size
#' r<-5 # number of variables
#' mis<-0.2 # frequency of missing data
#' mydata<-matrix(NA, nrow=n, ncol=r) # mydata is a matrix of r variables
#' # following a U(0,1) distribution
#' for (i in c(1:r)){
#' 	mydata[,i]<-runif(n,0,1)
#' }
#' bin.var<-sample(LETTERS[1:2],n,replace=TRUE, prob=c(0.3,0.7)) # binary variable [A,B].
#' # The probability of being in one of the categories is 0.3.
#' cat.var<-sample(LETTERS[1:3],n,replace=TRUE, prob=c(0.5,0.3,0.2)) # categorical variable [A,B,C].
# The vector of probabilities of occurence A, B and C is (0.5,0.3,0.7).
#' num.var<-runif(n,0,1) # Additional continuous variable following a U(0,1) distribution
#' mydata<-cbind.data.frame(mydata,bin.var,cat.var,num.var,stringsAsFactors = TRUE)
#' # dataframe with r+3 variables
#' colnames(mydata)=c("v1","v2","X1","X2","X3","X4","X5", "X6") # names of columns
#' # MCAR on X1 and X4 by using v1 and v2. MAR on X3 and X5 by using X2 and X6.
#' mydata$X1[which(mydata$v1<=sort(mydata$v1)[mis*n])]<-NA # X1: (mis*n)% of MCAR data.
#' # All data above the (100-mis)th percentile in v1 are selected
#' # and the corresponding observations in X1 are replaced with missing data.
#' mydata$X3[which(mydata$X2<=sort(mydata$X2)[mis*n])]<-NA # X3: (mis*n)% of MAR data.
#' # All data above the (100-mis)th percentile in X2 are selected
#' # and the corresponding observations in X3 are replaced with missing data.
#' mydata$X4[which(mydata$v2<=sort(mydata$v2)[mis*n])]<-NA # X4: (mis*n)% of MCAR data.
#' # All data above the (100-mis)th percentile in v2 are selected
#' # and the corresponding observations in X4 are replaced with missing data.
#' mydata$X5[which(mydata$X6<=sort(mydata$X6)[mis*n])]<-NA # X5: (mis*n)% of MAR data.
#' # All data above the (100-mis)th percentile in X6 are selected
#' # and the corresponding observations in X5 are replaced with missing data.
#' mydata$v1=NULL
#' mydata$v2=NULL

#' RBtest(mydata)
#'
#' @importFrom stats complete.cases ks.test lm p.adjust predict
#' @import nnet psych
#' @export

RBtest<-function(data){
	Type.1<-rep(NA,ncol(data))
	for (h in which(complete.cases(t(data))==FALSE)){ #only test the missing data on variables which contains missing data (logic)

		U=data.frame(data[,h],data[ ,colSums(is.na(data)) == 0],stringsAsFactors = TRUE) # we use for the first test only the complete variables

		#names(U1)[1] <- names(D)[h]

		new_A=data.frame(U[is.na(U[,1])==FALSE,],stringsAsFactors = TRUE) # "A" part (complete cases)
		new_B=data.frame(U[is.na(U[,1])==TRUE,],stringsAsFactors = TRUE)	# "B" part

		names(new_A)[1]="DV"

		# if the variable is binary (2 category) :

		if (is.factor(new_A$DV)==TRUE && length(levels(new_A$DV))==2) {
			reg_A=multinom(DV~.,trace=FALSE,data=new_A) #trace argument for silent computation.
			names(new_A)[1]=names(data)[h]
			u_obs_hat=reg_A$fitted.values
			u_obs_hat=data.frame(u_obs_hat, (1-u_obs_hat),stringsAsFactors = TRUE)
			u_mis_hat=predict(reg_A, new_B, "probs")
			u_mis_hat=data.frame(u_mis_hat, (1-u_mis_hat),stringsAsFactors = TRUE)

			#-- ks.test --#

			## BH adjustement!
			RES=p.adjust(c(ks.test(u_obs_hat[,1], u_mis_hat[,1])$p.value,
										 ks.test(u_obs_hat[,2], u_mis_hat[,2])$p.value),
									 method = "BH", n = 2)
			# >0.05 --> accept MCAR

			# And then assign the right value to Type.1

			if(RES[1]>0.05 && RES[2]>0.05) {
				Type.1[h]=0
			} else {
				Type.1[h]=1
			}
			names(U)[1] <- names(data)[h]
		}

		# if the variable is categorical with more than two categories :

		names(new_A)[1]="DV"
		if (is.factor(new_A$DV)==TRUE && length(levels(new_A$DV))>2) {
			reg_A=multinom(DV~.,trace=FALSE,data=new_A) #trace argument for silent computation.
			names(new_A)[1]=names(data)[h]
			u_obs_hat=reg_A$fitted.values
			u_mis_hat=predict(reg_A, new_B, "probs")

			#-- ks.test --#

			# First, compute the p-values of the ks.test for each category of the variable of interest
			KolSmir<-rep(NA,ncol(u_obs_hat))
			for (i in 1:ncol(u_obs_hat)){
				KolSmir[i]<-ks.test(u_obs_hat[,i], u_mis_hat[,i])$p.value
			}
			# And then --> BH adjustement!
			RES=p.adjust(KolSmir, method = "BH", n = ncol(u_obs_hat))
			# >0.05 --> accept MCAR
			# And then assign the right value to Type.1
			if (any(RES<0.05)==FALSE){
				Type.1[h]=0
			} else {
				Type.1[h]=1
			}
			names(U)[1] <- names(data)[h]
		}

		# if the variable is continuous :
		names(new_A)[1]="DV"
		if (is.factor(new_A$DV)==FALSE) {

			reg_A=lm(DV~.,data=new_A)
			names(new_A)[1]=names(data)[h]

			u_obs_hat=predict(reg_A, new_A) # prediction of x1_obs
			u_mis_hat=predict(reg_A, new_B) # prediction of x1_mis

			#-- ks.test --#
			RES=ks.test(u_obs_hat, u_mis_hat)$p.value
			# >0.05 --> accept MCAR

			if(RES>0.05) Type.1[h]=0
			if(RES<=0.05) Type.1[h]=1
			names(U)[1] <- names(data)[h]
		}
	}
	Type.1[which(is.na(Type.1))]=-1

	type<-Type.1
	names(type)<-names(data)
	abs_nbrMD<-colSums(is.na(data))
	rel_nbrMD<-colSums(is.na(data))/(nrow(data))
	return(list(abs.nbrMD=round(abs_nbrMD,0),rel.nbrMD=round(rel_nbrMD,2), type=type))

}

#' Test of missing data mechanism using all available information
#'
#' This function tests MCAR vs MAR by using both the complete and incomplete variables.
#' @usage RBtest.iter(data,K)
#' @param data Dataset with at least one complete variable.
#' @param K Maximum number of iterations.
#'
#' @return A list of the following elements:
#'  \itemize{
#'   \item \code{abs.nbrMD} The absolute quantity of missing data per variable.
#'   \item \code{rel.nbrMD} The percentage of missing data per variable.
#'   \item \code{K} The maximum admitted number of iterations.
#'   \item \code{iter} The final number of iterations.
#'   \item \code{type.final} Vector of the same length than the number of variables of the dataset, where '0' is for variables with MCAR data, '1' is for variables with MAR data and '-1' is for complete variables.
#'   \item \code{TYPE.k} Dataframe containing the type of missing data after each iteration. Each row is a vector having the same length than the number of variables of the dataset, where '0' is for variables with MCAR data, '1' is for variables with MAR data and '-1' is for complete variables.
#'    }
#'
#' @examples
#'
#' set.seed(60)
#' n<-100 # sample size
#' r<-5 # number of variables
#' mis<-0.2 # frequency of missing data
#' mydata<-matrix(NA, nrow=n, ncol=r) # mydata is a matrix of r variables
#' # following a U(0,1) distribution
#' for (i in c(1:r)){
#' 	mydata[,i]<-runif(n,0,1)
#' }
#' bin.var<-sample(LETTERS[1:2],n,replace=TRUE, prob=c(0.3,0.7)) # binary variable [A,B].
#' # The probability of being in one of the categories is 0.3.
#' cat.var<-sample(LETTERS[1:3],n,replace=TRUE, prob=c(0.5,0.3,0.2)) # categorical variable [A,B,C].
#' # The vector of probabilities of occurence A, B and C is (0.5,0.3,0.7).
#' num.var<-runif(n,0,1) # Additional continuous variable following a U(0,1) distribution
#' mydata<-cbind.data.frame(mydata,bin.var,cat.var,num.var,stringsAsFactors = TRUE)
#' # dataframe with r+3 variables
#' colnames(mydata)=c("v1","v2","X1","X2","X3","X4","X5", "X6") # names of columns
#' # MCAR on X1 and X4 by using v1 and v2. MAR on X3 and X5 by using X2 and X6.
#' mydata$X1[which(mydata$v1<=sort(mydata$v1)[mis*n])]<-NA # X1: (mis*n)% of MCAR data.
#' # All data above the (100-mis)th percentile in v1 are selected
#' # and the corresponding observations in X1 are replaced with missing data.
#' mydata$X3[which(mydata$X2<=sort(mydata$X2)[mis*n])]<-NA # X3: (mis*n)% of MAR data.
#' # All data above the (100-mis)th percentile in X2 are selected
#' # and the corresponding observations in X3 are replaced with missing data.
#' mydata$X4[which(mydata$v2<=sort(mydata$v2)[mis*n])]<-NA # X4: (mis*n)% of MCAR data.
#' # All data above the (100-mis)th percentile in v2 are selected
#' # and the corresponding observations in X4 are replaced with missing data.
#' mydata$X5[which(mydata$X6<=sort(mydata$X6)[mis*n])]<-NA # X5: (mis*n)% of MAR data.
#' # All data above the (100-mis)th percentile in X6 are selected
#' # and the corresponding observations in X5 are replaced with missing data.
#' mydata$v1=NULL
#' mydata$v2=NULL
#'
#' RBtest.iter(mydata,5)
#'
#' @importFrom stats complete.cases ks.test lm p.adjust predict
#' @import nnet mice psych
#' @export

RBtest.iter<-function(data,K){

	#Type.1 <- TestMDM(data)
	Type.1 <- RBtest(data)$type

	D1<-data
	TYPE_k=data.frame(matrix(NA, ncol=ncol(data)+1, nrow=K+1),stringsAsFactors = TRUE) #initialisation of the sequences of missing data mechanisms for each iteration

	Type.2<-rep(NA,ncol(data))
	k=1 # k is the number of iterations for the while loop. This is for excluding the case of infinite loop.

	TYPE_k[1,]=c(k,Type.1)
	Type.km1 <-rep(-2,(ncol(data)+1))

	while(all(TYPE_k[k,2:ncol(TYPE_k)]==Type.km1[2:length(Type.km1)])==FALSE) {

		for (h in which(complete.cases(t(data))==FALSE)){ #only test the missing data on variables which contains missing data (logic)

			U=data.frame(data[,h],D1[,-h],stringsAsFactors = TRUE) # first on v_h from the original dataset

			names(U)[1] <- names(data)[h]

			new_A=data.frame(U[is.na(U[,1])==FALSE,],stringsAsFactors = TRUE) # "A" part
			new_B=data.frame(U[is.na(U[,1])==TRUE,],stringsAsFactors = TRUE)	# "B" part

			names(new_A)[1]="DV"

			# if the variable is binary (2 category) :

			if (is.factor(new_A$DV)==TRUE && length(levels(new_A$DV))==2) {
				reg_A=multinom(DV~.,trace=FALSE,data=new_A) #trace argument for silent computation.
				names(new_A)[1]=names(data)[h]
				u_obs_hat=reg_A$fitted.values
				u_obs_hat=data.frame(u_obs_hat, (1-u_obs_hat),stringsAsFactors = TRUE)
				u_mis_hat=predict(reg_A, new_B, "probs")
				u_mis_hat=data.frame(u_mis_hat, (1-u_mis_hat),stringsAsFactors = TRUE)

				#-- ks.test --#

				## BH adjustement!
				RES=p.adjust(c(ks.test(u_obs_hat[,1], u_mis_hat[,1])$p.value,
											 ks.test(u_obs_hat[,2], u_mis_hat[,2])$p.value),
										 method = "BH", n = 2)
				# >0.05 --> accept MCAR

				# And then assign the right value to Type.1

				if(RES[1]>0.05 && RES[2]>0.05) {
					Type.2[h]=0
				} else {
					Type.2[h]=1
				}
				names(U)[1] <- names(data)[h]

				### Imputation of v_h with respect to the result!

				if (RES[1]>0.05 && RES[2]>0.05) {
					D1[is.na(data[,h]),h]<-sample(D1[!is.na(data[,h]),h],length(D1[is.na(data[,h]),h]),replace=TRUE)

				} else {
					U=complete(mice(U, m = 1, maxit=5, printFlag=FALSE))
					#printFlag option for silent computation
					D1[,h]=U[,1]
				}

			}
			names(new_A)[1]="DV"
			# if the variable is categorical with more than two categories :
			if (is.factor(new_A$DV)==TRUE && length(levels(new_A$DV))>2) {
				reg_A=multinom(DV~.,trace=FALSE,data=new_A) #trace argument for silent computation.
				names(new_A)[1]=names(data)[h]
				u_obs_hat=reg_A$fitted.values
				u_mis_hat=predict(reg_A, new_B, "probs")

				#-- ks.test --#

				# First, compute the p-values of the ks.test for each category of the variable of interest
				KolSmir<-rep(NA,ncol(u_obs_hat))
				for (i in 1:ncol(u_obs_hat)){
					KolSmir[i]<-ks.test(u_obs_hat[,i], u_mis_hat[,i])$p.value
				}
				# And then --> BH adjustement!
				RES=p.adjust(KolSmir, method = "BH", n = ncol(u_obs_hat))
				# >0.05 --> accept MCAR

				if (any(RES<0.05)==FALSE){
					Type.2[h]=0
				} else {
					Type.2[h]=1
					names(U)[1] <- names(data)[h]
				}

				# Imputation of v_h with respect to the result!

				if (any(RES<0.05)==FALSE) {
					D1[is.na(data[,h]),h]<-sample(D1[!is.na(data[,h]),h],length(D1[is.na(data[,h]),h]),replace=TRUE)

				} else {
					U=complete(mice(U, m = 1, maxit=5, printFlag=FALSE))
					#printFlag option for silent computation
					D1[,h]=U[,1]
				}
			}

			# if the variable is continuous :

			names(new_A)[1]="DV"

			if (is.factor(new_A$DV)==FALSE) {

				reg_A=lm(DV~.,data=new_A)
				names(new_A)[1]=names(data)[h]

				u_obs_hat=predict(reg_A, new_A) # prediction of x1_obs
				u_mis_hat=predict(reg_A, new_B) # prediction of x1_mis

				#-- ks.test --#
				RES=ks.test(u_obs_hat, u_mis_hat)$p.value
				# >0.05 --> accept MCAR

				if(RES>0.05) Type.2[h]=0
				if(RES<=0.05) Type.2[h]=1
				names(U)[1] <- names(data)[h]

				### Imputation of v_h with respect to the result!

				if (RES>0.05) {
					D1[is.na(data[,h]),h]<-sample(D1[!is.na(data[,h]),h],length(D1[is.na(data[,h]),h]),replace=TRUE)
				}
				if(RES<=0.05) {
					U=complete(mice(U, m = 1, printFlag=FALSE))
					#printFlag option for silent computation
					D1[,h]=U[,1]
				}
			}
		}


		Type.2[which(is.na(Type.2))]=-1

		k=k+1

		TYPE_k[k,]=c(k,Type.2)
		Type.km1 <- TYPE_k[k-1,]
		if (k==K+1) break
	}
	names(TYPE_k)[1]<-"k"

	abs_nbrMD<-colSums(is.na(data))
	rel_nbrMD<-colSums(is.na(data))/(nrow(data))
	TYPE.k<-TYPE_k[which(complete.cases(TYPE_k)==TRUE),]
	names(TYPE.k)<-c("k",names(data))
	type.final<-TYPE_k[nrow(TYPE.k),2:ncol(TYPE.k)]
	names(type.final)<-names(data)
	iter<-nrow(TYPE.k)
	return(list(abs.nbrMD=round(abs_nbrMD,0),rel.nbrMD=round(rel_nbrMD,2),K=K,iter=iter, type.final=type.final, TYPE.k=TYPE.k))
}
