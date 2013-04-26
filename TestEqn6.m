%PROGRAM TO SIMULATE DATA FOR THE PURPOSES OF ESTIMATING DIMSIONALITY.
%TAKEN FROM TestEqn4.m.
%
%PRETTY MUCH UNCHANGED EXCEPT (A) WE TYPICALLY LOOK AT FAR FEWER MUTATIONS
%AND (B) INSTEAD OF PLOTTING ANYTHING WE JUST BUILT THE MATRIX OF VALUES
%AND THEN CALL THE MLE CODE.

%COMMENTS IN CAPS MOSTLY BY DMW APR 2011. FIRST OBJECTIVE IS TO
%ONLY SAMPLE MUTATIONS ON THE PEAK; DO THIS BY SETTING s0 = 0 AND
%LEAVE b = 0. (THE LATTER SLIGHTLY A KLUDGE BUT OTHERWISE CODE ONLY WANTS
%MUTATIONS THAT MOVE TOWARD THE OPTIMUM FROM THE STARTING POINT
%WHICH AIN'T POSSIBLE IN LIGHT OF THE FORMER.)

%SECOND OBJECTIVE IS TO FLIP FIGURE 13 BY 90 DEGREES SO THAT 
%TRUE THETA IS ON X-AXIS AND ESTIMATE THETA IS ON Y-AXIS. THE IDEA
%IS THEN TO SUPERIMPOSE THE FRACTION OF IMAGINARY ESTIMATED
%THETAS ON THE SAME X-AXIS, AND MOST IMPORTANTLY, TO DISTINGUISH
%BETWEEN IMAGINARY THETAS WITH REAL PART ZERO AND THOSE WITH REAL
%PARTS 180.

%FINALLY ADD CODE TO VIOLATION ROTATION SYMMETRY AND MUTATIONAL ADDITIVITY
%IN PHENOTYPE SPACE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  DERIVED FROM APRIL 14 TestEqn4.m WITH JANUARY 2012 SIMULATOR
%  UPDATES MANUALLY FOLDED IN.
%
%  AND NOW, INSTEAD OF ACTUALLY BUILDING THE MATRIX OF COS(THETA) VALUES
%  AND COUNTING EIGENVALUES, LET'S JUST ACCUMULATE A SAMPLE OF COS(THETA)
%  VALUES AND TRY TO FIT IT TO THE EXPECTED DISTRIBUTION, AS PER REVIEWER
%  2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all;
% close all;

%CONSTANTS
s0 = 0.0; 
lambda_e = 10;  %This times next = s 
M = .1;%variance; stdev = sqroot of variance
%s IN THE PAPER IS lambda_e*M.
Selection = lambda_e*M;  %from Martin et al 2007; M = I_ne, I_ne is the identity matrix of dimension ne; it will be all 1's.


%input parameters

numdimensions = [5]; %vector with the #of dimensions
%ACTUALLY SINCE WE'RE NOW USING ALL MUTATIONS, THIS IS THE NUMBER OF
%MUTATION EXAMINED TOTAL.
numbenmutations = 25; %number of beneficial mutations to collect
Noise1 = 0.00;  %standard deviation in fitness measure
R = 1; %number of replicates for noisy fitness measures
a = 0.*Selection; %EXTENT OF NON-ROTATIONAL SYMMETRY IN SELECTION MATRIX S. 
PhenotypicEpistasis = 0*M; %Std deviation in phenotypic epistasis, to be added to each phenotype independently.
m11 = 0.0;
m12 = 0.0;



fignum = 1;

%ONCE UPON A TIME I PACKAGED THE WHOLE THING INSIDE A HUGE LOOP TO SAMPLE 
%THE VARIANCE IN THETAS

for replicate=1:1

%Flags
rnC = 0;%%this is the rounding # --MAKE SURE TO MAKE 0 IF YOU DON'T WANT TO ROUND
b = 0;  %this will allow you to switch between collecting all mutations (b = 0) or just beneficial mutations , b =1
SC = 0;  %this is the anticipated equal length for all the mutational vectors; keep at 0 to not force equal lengths
MAP = 0; %this flag will maps the relationship between scaled episatsis for all pairs of mutations - it is computationally intensive so make 0-

%%make arrays that data will be collected into
propben = ones(1, length(numdimensions));
sben = ones(length(numbenmutations),length(numdimensions));
theta = ones(length(numbenmutations),length(numdimensions));
r = ones(length(numbenmutations),length(numdimensions));
mean_Epistasis = ones(length(numdimensions),1);
var_Epistasis =ones(length(numdimensions),1);
beneficialtheta_compensatory = [];   ScaledEpistasis_compensatory = []; ratiolength_compensatory = [];
beneficialtheta_Ncompensatory = []; ScaledEpistasis_Ncompensatory = []; ratiolength_Ncompensatory = [];
beneficialtheta_Pcompensatory = []; ScaledEpistasis_Pcompensatory = []; ratiolength_Pcompensatory = [];
beneficialtheta_lesszero = []; ScaledEpistasis_lesszero=[];

%determine the number of Epistasis measurements and make an array to
%collect these data
C = (nchoosek(numbenmutations,2));
Epistasis = ones(C,size(numdimensions,2));
EpistasisPredictsTheta = ones(C,size(numdimensions,2)); 
PredictTheta = ones(C,size(numdimensions,2));   
ScaledEpistasis = ones(C,size(numdimensions,2));
beneficialtheta = ones(C,size(numdimensions,2));
ratiolength = ones(C,size(numdimensions,2));
CompareTheta_RatioChange = ones(C,3);
CompareTheta = ones(C,2);
sumbendzi = ones(C,numdimensions);

%make some vectors to collect noisy fitness measures- troubleshooting
WxNoisyAll = ones(numbenmutations,length(numdimensions));
WyNoisyAll =ones(numbenmutations,length(numdimensions));
WxyNoisyAll = ones(numbenmutations,length(numdimensions));
WxAll = ones(numbenmutations,length(numdimensions));
WyAll = ones(numbenmutations,length(numdimensions));
WxyAll = ones(numbenmutations,length(numdimensions));
WxNoisyProblem = []; WyNoisyProblem = []; WxyNoisyProblem = [];
CompareWxProblem = []; CompareWyProblem = []; CompareWxyProblem = [];


%This loop goes through each element in the dimension vector and for each
%n, it finds the given number of beneficial mutations, and collects
%information we care about, and then calculates Epistasis among all pairs
%of the beneficial mutations
for z = 1: length(numdimensions)
    n = numdimensions(z);

Selection_array = repmat(Selection,1,n);  
M_array = repmat(M,1,n);
sigma = diag(M_array); %M = covariance matrix


if (m11 ~= 0.)
    sigma(1,1) = m11;
end
if (m12 ~= 0.)
    sigma(1,2) = m12;
    sigma(2,1) = m12;
end

S = diag(Selection_array);
%IMPLEMENT NON-ROTATIONAL SYMMETRY IN SELECTION MATRIX
%if (n == 2) 
%    S(2,2) = S(2,2)*a;
%end

% JUNE 2012: THIS IS A BIT KLUDGY: NOT EVERY SYMMETRIC MATRIX IS
% SEMI-DEFINITE SO WE'LL JUST KEEP TRYING TILL WE FIND ONE!

stemp = S;
while 1 == 1
    S = stemp + randn(n)*a;

    for i=1:n
        for j=i+1:n
            S(j,i) = S(i,j);
        end
    end
    if min(real(eig(S))) >= 0 break
    end
end

%PLAYING AROUND WITH NON-ZERO NON-MAIN DIAGONAL
%S(1,2) = S(1,2)*1.5;
%S(2,1) = S(2,1)*2;
%S(2,2) = S(2,2)*3;

%% find the individual trait values of the starting genotype and calculate
%% the fitness of this genotype
z0_individual = (2*s0/(n*Selection))^0.5;
    if rnC > 0
    z0_individual = round(z0_individual*rnC)/rnC;
    else
    end
z0 = repmat(z0_individual,n,1);             %intitial distance to the optimum scaled to equal s0
W_z0 = exp(-0.5*z0'*S*z0) ;                  %fitness of initial genotype 

totalmutations = 0;                         % initialize the mutation counter
m = 1;                                      % initialize the beneficial mutation counter
bendzi = ones(numbenmutations,n);
bendziNoisy = ones(numbenmutations,n);

%this outer loop is collecting random mutations - the inner if loop is
%collecting interesting data for only the beneficial mutations
while m <= numbenmutations  
%   Jen used mgd() here but chatter on the web makes me nervious about
%   using it with non-zero covariance. Calling randn() directly seems to be
%   the recommended approach: http://www.mathworks.com/help/matlab/ref/randn.html
%
%     x = mgd(1,n,z0,sigma);               %random draw of numbenmutations from multivariate Gaussian distribution of r z0 vectors

    x = randn(1,n);
    x = x*sigma;
    
    if rnC >0
    x = round(x.*rnC)./rnC;
    else 
    end
    
    x = x' ;                              % transform into a column vector
    dzi = x - z0;                           %calculate dzi - the change in each phenotype per mutation x
    
    if SC ==0
    else
     c = SC/(sum(dzi.^2)^0.5);              %this is scaling the mutations to all have the same length
     dzi = dzi.*c;         
    end      
       
    W_z_random= exp(-0.5*x'*S*x);           %fitness of the mutant
    log_W_z_random = log(W_z_random/W_z0);  %this is one way to calculate the log fitness of the mutant
    log_wi = -Selection*(z0'*dzi+0.5*dzi'*dzi); %and this another. this is actually a different way of expressing the first. you get the same answer with either.
    s = (W_z_random./W_z0)-1;                 %calculate the selection coefficient of the mutation
        
        
        %for beneficial mutations (those with s >0), this loop will: 1)
        %store the dzi, 2) store their selection coefficient, 3) find the
        %angle of the mutation relative to the angle going straight to the
        %optimum and store this angle
        
        
        if (b == 1) 
            if log(1+s) <= 0
            else   %this is a flag to just collect beneficial mutations
            bendzi(m,:) = dzi';
            sben(m,z) = s;  % keep the 's' for each beneficial mutation
            %find the angle - theta - of each beneficial mutation- refer to
            %Poon and Otto 2001, equation 19, solving for theta; x =
            %cos(theta); x here is (z'^2- z^2-r^2)/(-2zr) 
            newz_2 = sum(x.^2);
            z_2 = sum(z0.^2);
            originalz = z_2^0.5;
            r_2 = sum(dzi.^2);
            r(m,z) = r_2^0.5;
                
            
            angleinput = (newz_2-z_2-r_2)/(-2*originalz*r(m,z));
            theta(m,z) = acosd(angleinput);
            m = m+1;
            end
        elseif b ==0 
            bendzi(m,:) = dzi';
            sben(m,z) = s;  % keep the 's' for each beneficial mutation
            %find the angle - theta - of each beneficial mutation- refer to
            %Poon and Otto 2001, equation 19, solving for theta; x =
            %cos(theta); x here is (z'^2- z^2-r^2)/(-2zr) 
            newz_2 = sum(x.^2);
            z_2 = sum(z0.^2);
            originalz = z_2^0.5;
            r_2 = sum(dzi.^2);
            r(m,z) = r_2^0.5;
            angleinput = (newz_2-z_2-r_2)/(-2*originalz*r(m,z));
            theta(m,z) = acosd(angleinput);
            m = m+1;
        else
        end

        %keep track of the proportion of beneficial mutations         
        totalmutations = totalmutations+1;
end
    
propben(1,z) = numbenmutations/totalmutations;
   

%calculate Epistasis and the angle between all possible pairs of beneficial mutations ON THE
%WILDTYPE BACKGROUND - where the mutations are now deleterious!
EpistasisNoisy = []; PredictThetaNoisy = [];  EpistasisPredictsThetaNoisy =[]; PredictThetaCompare =[]; EpistasisNoisyFiltered=[];
SEpistasis_Array = tril(ones(numbenmutations,numbenmutations),0);
Epistasis_Array = tril(ones(numbenmutations,numbenmutations),0);
greater0 = 0;
lesszero = 0;
greater0P = 0;
greater0N = 0;
l =1; j = 2; k = 1;
CountImagInvTheta = 0;

%KLUDGE-TIME! I WANT A MATRIX TO HOLD ONLY THE REAL ESTIAMTED THETAS
%TOGETHER WITH THEIR UNDERLYING TRUE THETAS SO'S I CAN COMPUTE r^2.
%KLUDGE BECAUSE I'M JUST BOLTING THIS ON TOP OF EXISTING CODE INSTEAD OF
%THINKING ABOUT AN ELEGANT WAY TO PULL JUST THESE VALUES AT THE BACK END.
%MORE KLUDGE -- ALTHOUTH ALL FIGURES DRAWN OVER 200 MUTATIONS, r^2 COMPUTES
%OVER ONLY 100 B/C regstats() BLEW UP WHEN MATRICIES HAD 50*99 ENTRIES.
%UBER-KLUDGE -- WILL WRITE THIS TO A FILE FOR EXCEL TO READ BECAUSE I CAN'T
%FIGURE OUT HOW TO FORCE THE REGRESSION THROUGH THE ORIGIN...!
%
%TO MAKE THIS WORK RUN THIS FILE AND THEN FROM THE MATLAB WINDOW TYPE
%
%save('rsst','rSquare','-ASCII') [HEY DUMMY: OMIT THAT LEADING '%'!)

%THEN OPEN rsst (WHICH RESIDES IN THE USER'S MATLAB FOLDER) WITH EXCEL,
%SPECIFYING SPACE AS THE DELIMITER. SORT THE TWO COLUMNS TOGETHER AND
%DELETE THE ZEROS (NOTE THAT WE INITIALIZE THE ARRAY WITH zero()! SO THESE 
%ARE WHAT'S LEFT OVER AS A CONSEQUENCE OF IMAGINARY THETAS), %SCATTERPLOT, 
%DRAW A LINEAR TRENDLINE THROUGH THE ORIGIN AND GRAB r^2.
%
%rSquareTrueTheta = [];
%rSquareEstTheta = [];
rSquare = zeros(C,2);

for l = 1:C
%     Epistasis(l,z) = -1*bendzi(k,:)*Selection*bendzi(j,:)'; %Epistasis is calculated according to equation 2 in Martin et al 2007
    
       %find the angle between all pairs of beneficial mutations, 1) find the
    %length of the vector between the 2 beneficial mutations
    %lenbenside = (sum((bendzi(k,:)-bendzi(j,:)).^2))^0.5;
    %beneficialtheta(l,z) = acosd((lenbenside^2-(r(k,z))^2-(r(j,z))^2)./(-2*r(k,z)*r(j,z)));
% The above code works if we don't violate model assumptions. But the more
% robust way is to use the dot product.

% THIS IS THE WAY WE HAD IT CODED FOR THE ORIGINAL SUBMISSION (WITH THE
% EXCEPTION THAT WE USED VALUES PREVIOUSLY STORED IN r(), WHICH IS ONLY
% KOSHER IF z0 = 0).
%
%    beneficialtheta(l,z) = acosd(dot(bendzi(j,:),bendzi(k,:))/sqrt(dot(bendzi(j,:),bendzi(j,:))*dot(bendzi(k,:),bendzi(k,:))));

% JANUARY 2, 2012: AS PER REVIEWERS' SUGGESTION, THIS MODEL MORE ROBUST
% THAN WE UNDERSTOOD. SHOULD WORK FINE FOR ARBITRARY GAUSSIAN FITNESS
% FUNCTION, EXCEPT THAT THE *TRUE* UNDERLYING THETA MUST SIMILARLY BE
% COMPUTED USING TRUE MATRIX S VIA THE INNER PRODUCTS x*S*y and x*S*x 
% INSTEAD OF DOT PRODUCTS.
%
% NOTE THAT THE ORIENTATION OF THE bendzi VECTORS SEEMS TO BE BACKWARDS
% (THEY'RE ROWS, NOT COLUMNS) WITH THE CONSEQUENCE THAT I HAVE TO TRANSPOSE
% THE LATTER RATHER THAN THE FORMER IN EACH INNER PRODUCT.

    beneficialtheta(l,z) = acosd(bendzi(j,:)*(S*bendzi(k,:)')/(sqrt(bendzi(j,:)*(S*bendzi(j,:)'))*sqrt(bendzi(k,:)*(S*bendzi(k,:)'))));
 
    CompareTheta_RatioChange (l,:) = [beneficialtheta(l,z) bendzi(j,1)./bendzi(j,2) bendzi(k,1)./bendzi(k,2)];
    CompareTheta(l,:) = [theta(j) theta(k)];
    
  
     
        worstcase = zeros(z,1);
        worstcase(1,1) = r(j,z)+r(k,z);
             


        fitnessWxy = fitness(S,z0+bendzi(j,:)'+bendzi(k,:)'+randn(n,1)*PhenotypicEpistasis);
        WxyAll(l,z) = fitnessWxy;
        fitnessWx = fitness(S,z0+bendzi(j,:)');
             WxAll(l,z) = fitnessWx;
        fitnessWy = fitness(S,z0+bendzi(k,:)');
             WyAll(l,z) = fitnessWy;
        fitnessWxWy = fitness(S,bendzi(j,:)')*fitness(S,bendzi(k,:)');
        fitnessWxWyP = fitness(S,worstcase);
        fitnessWxWyP_Dan = exp(-0.5*(-2*log(fitnessWx)-2*log(fitnessWy)+4*(log(fitnessWx)*log(fitnessWy))^0.5));
        
        
     Epistasis(l,z) = fitness(S,bendzi(j,:)'+bendzi(k,:)')-(fitness(S,bendzi(j,:)')*fitness(S,bendzi(k,:)'));
     Epistasis_Array(k,j) = Epistasis(l,z);
     test1 = log((Epistasis(l,z)/(exp(log(fitness(S,bendzi(j,:)')))*exp(log(fitness(S,bendzi(k,:)')))))+1);
%     EpistasisPredictsTheta(l,z) = log((Epistasis(l,z)/(exp(log(fitness(S,bendzi(j,:)')))*exp(log(fitness(S,bendzi(k,:)')))))+1)/(-2*(((log(fitness(S,bendzi(j,:)')))*(log(fitness(S,bendzi(k,:)'))))^0.5));
     EpistasisPredictsTheta(l,z) = Eqn3(fitnessWx,fitnessWy,fitnessWxy);
%     PredictTheta(l,z) = acosd(EpistasisPredictsTheta(l,z));
% ABOVE ONLY WORKS IF WE DON'T VIOLATE MODEL ASSUMPTIONS. BELOW IS MORE
% ROBUST: USES THE DOT PRODUCT.
    PredictTheta(l,z) = beneficialtheta(l,z);

%      WxyNoisy = mgd(1,1,fitnessWxy,Noise);
%      WyNoisy = mgd(1,1,fitnessWy,Noise);  
%      WxNoisy =mgd(1,1,fitnessWx,Noise);
     
     
     WxyNoisy = mean(randn(R,1).*Noise1+fitnessWxy); WxyNoisyAll(l,z) = WxyNoisy;
     WyNoisy = mean(randn(R,1).*Noise1+fitnessWy);  WyNoisyAll(l,z) = WyNoisy;
     WxNoisy =mean(randn(R,1).*Noise1+fitnessWx);    WxNoisyAll(l,z) = WxNoisy;
     
 % JAN 2, 2012: NEED THIS NOW THAT WE ALLOW BACKGROUND PHENOTYPE TO BE 
 % SOMEWHERE OTHER THAN THE ORIGIN.
      if (s0 > 0) W0Noisy = mean(randn(R,1).*Noise1+fitness(S,z0));
      else W0Noisy = 1.;
      end     

     % THIS CODE PROTECTS US FROM CASES IN WHICH THE LOG RETURNS AN
     % IMAGINARY NUMBER, IN WHICH CASE WE CAN'T EVEN TAKE acosd().
 
     % BUT IN FACT WHAT WE DO WITH THE REAL DATA IS TO NORMALIZE SO THAT
     % THE LARGEST (NOISY) FITNESS VALUE IS FIXED AT 1. SO WHY DON'T WE DO
     % THAT HERE? AFTER ALL THE ONLY PLACE I CAN SEE AN IMAGINARY NUMBER
     % CREEPING IN IS IF WE TAKE THE sqrt() OF THE PRODUCT OF TWO LOGS OF
     % OPPOSITE SIGN, I.E IF ONE (BUT NOT BOTH) OF THE NOISY SINGLE MUTANTS
     % IS GREATER THAN 1.0
     
     % OOPS: ONE MORE PLACE AN IMAGINARY NUMBER CAN APPEAR: IF AFTER ADDING
     % NOISE ANY MUTANTS HAVE FITNESS OF ZERO OR WORSE.
     
     % SLIGHTLY KLUDGY FIX: JUST SHAVE THEM DOWN TO BE EPSILON LESS THAN 1
     % OR SHAVE THEM UP TO BE EPSILON GREATER THAN ZERO.
     % IN THE REAL DATA WE WOULDN'T EVEN CONSIDER SUCH A PAIR SINCE THEY'RE
     % NOT BOTH DELETERIOUS. SO IF WE HAD IT TO DO OVER AGAIN HERE WE MIGHT
     % REDRAW UNDERLYING MUTATIONS, BUT I DON'T WANT TO BOTHER WITH THAT
     % NOW.

     if (WxNoisy > 1) WxNoisy = 0.999999999999999;CountImagInvTheta=CountImagInvTheta+1;
     end
     if (WyNoisy > 1) WyNoisy = 0.999999999999999;CountImagInvTheta=CountImagInvTheta+1;
     end
     if (WxNoisy < 0) WxNoisy = 0.000000000000001;CountImagInvTheta=CountImagInvTheta+1;
     end
     if (WyNoisy < 0) WyNoisy = 0.000000000000001;CountImagInvTheta=CountImagInvTheta+1;
     end
     if (WxyNoisy < 0) WxyNoisy = 0.000000000000001;CountImagInvTheta=CountImagInvTheta+1;
     end
     if (WxyNoisy> 1) 
     end
     
     WxWyPNoisy = exp(-0.5*(-2*log(WxNoisy)-2*log(WyNoisy)+4*(log(WxNoisy)*log(WyNoisy))^0.5));
%     if (WxNoisy <1 && WyNoisy <1 && WxyNoisy <1 && WxyNoisy > WxWyPNoisy && WxNoisy > WxWyPNoisy && WyNoisy > WxWyPNoisy)
        EpistasisNoisy= [EpistasisNoisy (WxyNoisy-(WxNoisy*WyNoisy))];
        %JEN'S EQUATION HERE: WHAT'S GOING ON WITH
        %THE +1 IN THE NUMERATOR? (THERE'S ALSO A COUPLE exp(log(x))'S BUT
        %I DON'T WORRY AS MUCH ABOUT THOSE.
        NoisyInverseTheta = log(WxyNoisy*W0Noisy/(WxNoisy*WyNoisy))/(-2*((log(WxNoisy)*log(WyNoisy))^0.5));
 %        NoisyInverseTheta =  log((WxyNoisy-(WxNoisy*WyNoisy))/(exp(log(WxNoisy))*exp(log(WyNoisy)))+1)/(-2*(((log(WxNoisy))*(log(WyNoisy)))^0.5));
 %    NoisyInverseTheta = Eqn3(WxNoisy,WyNoisy,WxyNoisy);
            if (abs(NoisyInverseTheta) >1)
                WxNoisyProblem = [WxNoisyProblem WxNoisy];
                WyNoisyProblem = [WyNoisyProblem WyNoisy];
                WxyNoisyProblem = [WxyNoisyProblem WxyNoisy];
           
                CompareWxProblem = [CompareWxProblem fitnessWx];
                CompareWyProblem = [CompareWyProblem fitnessWy];
                CompareWxyProblem = [CompareWxyProblem fitnessWxy];
                
            else
                EpistasisNoisyFiltered = [EpistasisNoisyFiltered (WxyNoisy-(WxNoisy*WyNoisy))];
                EpistasisPredictsThetaNoisy = [EpistasisPredictsThetaNoisy NoisyInverseTheta];
            end
            
 % THE ABOVE else CLAUSE USED TO INCLUDE CODE TO LOAD PredictThetaCompare 
 % AND PredictThetaNoisy.
 % BUT NOW I WANT TO ALSO RECORD IMAGINARY THETAS IN PredictThetaNoisy[]
 % SINCE WE PROTECT OURSELVES, acosd() WILL ONLY EVER RETURN VALUES ON 
 % [0,180] SO I CAN CODE THIS AS -1 FOR 0+xi AND 181 FOR 180+xi. 
 
            PredictThetaCompare = [PredictThetaCompare PredictTheta(l,z)];
            
 % JULY 3, 2012: WE WANT TO DROP ALL IMAGINARY THETAS

            if ((NoisyInverseTheta < 1) & (NoisyInverseTheta > -1))
                PredictThetaNoisy = [PredictThetaNoisy acosd(NoisyInverseTheta)];
            end
            
        

    j = j+1;
    if j > numbenmutations
        k = k+1;
        j = k+1;
    end
end
  
% JUNE 28, 2012 AT THIS POINT PredictedNoisyTheta() HAS OUR ESTIMATED
% THETAS. NOW LET'S TRY TO FIND WHAT n FITS THESE DATA BEST.

% ThetaLogLikelihood() wants its angles in radians, so convert them here.

PredictThetaNoisy = PredictThetaNoisy * pi / 180.;

x = zeros(C,1);
p = zeros(C,1);

MAXp = -10000000000000000;
SUMp = 0;
MINn = 1;
MAXn = 40;
nSTEP = 0.05;
x = zeros((MAXn-MINn)/nSTEP,1);
p = zeros((MAXn-MINn)/nSTEP,1);
for i=1:(MAXn-MINn)/nSTEP
    x(i) = (i-1)*nSTEP + MINn;
    % Matrix needs to be transposed
    p(i) = ThetaLikelihood(x(i),PredictThetaNoisy');
    SUMp = SUMp + p(i);
    if p(i) > MAXp
        MAXp = p(i);
        CI(1) = x(i);
    end
end

lowflag = 0;
highflag = 0;
runningSUMp = 0;
for i=1:(MAXn-MINn)/nSTEP
    runningSUMp = runningSUMp + p(i);
    if ((runningSUMp/SUMp > 0.025) & (lowflag == 0))
        CI(2) = x(i);
        lowflag = 1;
    end
    if ((runningSUMp/SUMp > 0.975) & (highflag == 0))
        CI(3) = x(i);
        highflag = 1;
    end
end

plot(x,p)
CI(4) = size(PredictThetaNoisy,2);
CI

n = CI(1);
for(i=1:181)
    null_cdf(i,1) = i-1;
    s = sincdf(i-1,n);
    % Annoying behavior inside null_cdf() that behaves badly very near 0.0 
    % and very near 1.0!
    if (s < 1e-12)
        null_cdf(i,2) = 0.;
    elseif (abs(1 - s) < 1e-12)
        null_cdf(i,2) = 1.;
    else
        null_cdf(i,2) = s;
    end
end
[h0,p,ksstat,cv]=kstest(PredictThetaNoisy*180/pi,null_cdf,0.05)


end


end
