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

clear all;
% close all;

%constants
s0 = 0.0; 
lambda_e = 10;  %This times next = s 
M = .1;%variance; stdev = sqroot of variance
%s IN THE PAPER IS lambda_e*M.
Selection = lambda_e*M;  %from Martin et al 2007; M = I_ne, I_ne is the identity matrix of dimension ne; it will be all 1's.

%SET INPUT PARAMETERS HERE! ##########################################
%
numdimensions = [20]; %vector with the #of dimensions. Though not sure it works with more than one element...
numbenmutations = 50; %number of beneficial mutations to collect
Noise1 = 0.00;  %standard deviation in fitness measure
Noise2 = 0;
R = 1; %number of replicates for noisy fitness measures
a = 0; %ADD A NORMALLY DISTR'D MATRIX WITH SD a TO SELECTION MATRIX S. This never matters, as our algebra is robust to deviations from identity matrix
PhenotypicEpistasis = 0.0; %Std deviation in phenotypic epistasis, to be added to each phenotype independently.
%
%END OF INPUT PARAMETERS #############################################


%Flags
rnC = 0;%%this is the rounding # --MAKE SURE TO MAKE 0 IF YOU DON'T WANT TO ROUND
b = 0;  %this will allow you to switch between collecting all mutations (b = 0) or just beneficial mutations , b =1
SC = 0;  %this is the anticipated equal length for all the mutational vectors; keep at 0 to not force equal lengths
MAP = 0; %this flag will maps the relationship between scaled episatsis for all pairs of mutations - it is computationally intensive so make 0-

%%make arrays that data will be collected into. 
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
% sigma(1,2) = 0.099;
S = diag(Selection_array);
%IMPLEMENT NON-ROTATIONAL SYMMETRY

S = S + randn(n)*a;

% BUT MATRIX MUST ITSELF REMAIN SYMMETRIC!
% THIS WAS COMMENTED OUT IN APRIL VERSION; UNCOMMENTED JAN 2, 2012

for i=1:n
    for j=i+1:n
        S(j,i) = S(i,j);
    end
end

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

xSx = zeros(numbenmutations,2);

%this outer loop is collecting random mutations - the inner if loop is
%collecting interesting data for only the beneficial mutations
while m <= numbenmutations    
    x = mgd(1,n,z0,sigma);               % random draw of numbenmutations from 
                                         % multivariate Gaussian distribution of r z0 vectors
    
    
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
            if log(1+s) >= 0
            else % this means we only collect deleterious mutations.
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
            
            xSx(m,1) = sqrt(-2*log(fitness(S,z0+bendzi(m,:)')));
            xSx(m,2) = sqrt(bendzi(m,:)*S*bendzi(m,:)');
            
            m = m+1;
            end
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
%DELETE THE ZEROS (NOTE THAT WE INITIALIZE THE ARRAY WITH zero()!),
%SCATTERPLOT, DRAW A LINEAR TRENDLINE THROUGH THE ORIGIN AND GRAB r^2.
%
%rSquareTrueTheta = [];
%rSquareEstTheta = [];
rSquare = zeros(C,2);

% I WILL USE THIS TO STORE BOTH sqrt(-2*Log(Wx)) AND xSx. WHEN s0 = 0 THEY
% SHOULD BE IDENTICAL. OTHERWISE WE CAN SEE HOW BADLY THEY DIVERGE.

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
             
% JAN 2, 2012: WHEN I TESTED THE INNER PRODUCT CODE WITH NON-SYMMETRIC S,
% IT DIDN'T IMMEDIATELY WORK (BUMMER!) BUT THE ALGEBRA IS NOW SO SIMPLE
% THAT IT WAS EASY TO DEBUG. NAMELY, I MANUALLY COMPUTED EACH OF THE INNER
% PRODUCTS, AND THEN MANUALLY COMPUTED EACH OF THE CORRESPONDING ELEMENTS
% IN EQUATION 4. 
%
% TURNED OUT THE DISCREPANCY WAS BETWEEN INNER-PRODUCT-COMPUTED NORM OF
% bendzi(k,:) AND fitnessWy, WHICH IDENTIFIED A BUG IN COMPUTING THE
% LATTER. NAMELY, CODE USED "Selection" INSTEAD OF "S" AS THE ARGUMENT TO
% fitness().
%
% SAME BUG PROPOGATED A FEW MORE TIMES... NOT CLEAR WHAT EFFECT THIS HAD ON
% PREVIOUS SIMULATION RESULTS THOUGH IT SEEMS TO ONLY MATTER WHEN SELECTION 
% MATRIX IS NOT SYMMETRIC.

% Jan 4, 2012: AND HERE'S SOMETHING UNRELATED BUT TREMENDOUSLY IMPORTANT:
% WHEN WE'RE NOT PUTTING MUTATIONS ON THE ORIGIN, WE NEED TO ADD z0 TO EACH
% dz BEFORE COMPUTING FITNESS.

        fitnessWxy = fitness(S,z0+bendzi(j,:)'+bendzi(k,:)'+randn(n,1)*PhenotypicEpistasis);
        WxyAll(l,z) = fitnessWxy;
        fitnessWx = fitness(S,z0+bendzi(j,:)');
             WxAll(l,z) = fitnessWx;
%       fitnessWy = fitness(Selection,bendzi(k,:)');
        fitnessWy = fitness(S,z0+bendzi(k,:)');

             WyAll(l,z) = fitnessWy;
        fitnessWxWy = fitness(S,bendzi(j,:)')*fitness(S,bendzi(k,:)');
%       fitnessWxWyP = fitness(Selection,worstcase);
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
%      WxyNoisy = WxyNoisy / W0Noisy;
%      WxNoisy = WxNoisy / W0Noisy;
     
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
            if (NoisyInverseTheta > 1)
                PredictThetaNoisy = [PredictThetaNoisy -1];
            else if (NoisyInverseTheta < -1)
                    PredictThetaNoisy = [PredictThetaNoisy 181];
                else
                    PredictThetaNoisy = [PredictThetaNoisy acosd(NoisyInverseTheta)];
                    rSquare(l,1) = PredictTheta(l,z);
                    rSquare(l,2) = acosd(NoisyInverseTheta);
                end
            end
            

        

    j = j+1;
    if j > numbenmutations
        k = k+1;
        j = k+1;
    end
end
   
end

%stats = regstats(rSquareEstTheta,rSquareTrueTheta);

%data we care about
mean_sben = mean(sben);
mean_Epistasis = mean(Epistasis)';
var_Epistasis = var(Epistasis)';

if MAP ==1
[row,col] = phenotype(SEpistasis_Array);
else
end


    % HERE'S FIGURE 13 ROTATED 90 AND WITH FREQUENCIES OF BOTH KINDS OF
    % IMAGINARY THETAS OVERLAYED.
    
         num = length(PredictThetaCompare);
     
     
     % compute and plot sliding run mean +/- 1 stdev
    Thetas = zeros(num,2);
    SortedThetas = zeros(num,2);
    MeanVar = zeros(num,2);
    HighNums = zeros(num,1);
    Imaginaries = zeros(num,2);
    %I'LL STORE THE TWO KINDS OF IMAGINARY COUNTS IN EACH BIN HERE, ALONG
    %WITH THE NUMBER OF OBSERVATIONS TOTAL IN THE BIN. AND I AS I ENCOUNTER
    %IMAGINARY THETAS I HAVE TO COUNT THEM (BOTH TYPES). TO GET mean() AND 
    %std() TO WORK CORRECTLY I COPY THE OBSERVED REALS INTO A DUMMY ARRAY.
    %I CAN'T SIMPLY SHIFT-LEFT WITHIN Thetas() BECAUSE I NEED TO KEEP THE
    %IMAGINARIES AROUND TO BE COUNTED IN THE NEXT BIN.
    
    Thetas(1:num,1) = PredictThetaNoisy;
    Thetas(1:num,2) = PredictThetaCompare;
    % SORT BY COLUMN 2 INSTEAD OF 1
    SortedThetas = sortrows(Thetas,2);
 %   SortedThetas(1:10,:)
    for index=1:num
        % the interval between index and high captures 2.5 degrees of
        % realized (true) thetas. EXCEPT THAT I'VE DECIDED I LIKE A SLIDING
        % WINDOW 10 DEGREES WIDE!
        high = index;
        while (SortedThetas(high,2) < SortedThetas(index,2) + 10), 
            high = high+1;
            if (high > num) 
                high = num;
                break;
            end
        end
        % FOR EACH POSITION ON THE X-AXIS WE HAVE TO COUNT THE NUMBER OF
        % EACH OF THE TWO KINDS OF IMAGINARIES, AND COPY THE NON-IMAGINARIES
        % INTO A DUMMY HOLDER FOR mean() AND std() TO PROCESS. CAN'T BE ANY
        % BIGGER THAN high-index, WHICH IS THE WIDTH OF THIS BIN. BUT WE
        % NEED ONE MORE INDEX, SINCE THERE MAY BE FEWER THAN THAT NUMBER
        % OF REALS.
        DummyThetaHolder = zeros(high-index);
        DTHindex=0;
        for subindex=index:high
            if ((SortedThetas(subindex,1) == -1) || (SortedThetas(subindex,1) == 181))
                if (SortedThetas(subindex,1) == -1) Imaginaries(index,1)=Imaginaries(index,1)+1;
                else Imaginaries(index,2)=Imaginaries(index,2)+1;
                end
            else 
                DTHindex = DTHindex + 1;
                DummyThetaHolder(DTHindex) = SortedThetas(subindex,1);
            end
        end
        if (DTHindex ~= high-index+1)
        end
        MeanVar(index,1) = mean(DummyThetaHolder(1:DTHindex));
        MeanVar(index,2) = std(DummyThetaHolder(1:DTHindex));
        %DIVIDE OBSERVED IMAGINARIES BY TOTAL NUMBER OF OBSERVATIONS.
        Imaginaries(index,1) = Imaginaries(index,1)/(high-index);
        Imaginaries(index,2) = Imaginaries(index,2)/(high-index);
        % don't have any idea what this was once used for. Now never
        % referenced.
        HighNums(index) = high;
    end
    figure (7)
% SUBTRACT 5 FROM THE X-AXIS TO MAKE UP FOR THE FACT THAT THE BINS LIE 10
% DEGREES TO THE RIGHT.
%
% FINALLY, SCALE THIS IMAGE BY 55% TO FIT INTO FINAL FIGURE
plot(SortedThetas(:,2),MeanVar(:,1)+MeanVar(:,2)-5,'k-','LineWidth', 1.5);
hold on;
%set(gca,'DataAspectRatioMode','manual')
set(gca,'DataAspectRatio',[1,1,1])
set(gca,'FontSize',16)
set(gca,'XTick',0:20:180)
%set(gca,'XTick',0:20:180)
set(gca,'XLim',[0,180])
set(gca,'YTick',0:20:180)
set(gca,'YLim',[0,180])
%axis square
plot([0,180],[0,180],'k--','LineWidth',2)
[haxes,hline1,hline2]=plotyy(SortedThetas(:,2),MeanVar(:,1)-5,SortedThetas(:,2),Imaginaries(:,2));
set(haxes(1),'FontSize',16)
set(haxes(2),'FontSize',16)
set(haxes(1),'XTick',0:20:180)
set(haxes(2),'XTick',0:20:180)
set(haxes(1),'YLim',[0,180])
set(haxes(1),'YTick',0:20:180)
set(haxes(2),'YLim',[0,1])
set(haxes(2),'YTick',0:.2:1)
set(hline1,'LineStyle','-')
set(hline1,'LineWidth',1.5)
set(hline1,'Color','k')
set(haxes(1),'YColor','k')
set(hline2,'LineStyle','--')
set(hline2,'LineWidth',.5)
set(hline2,'Color','k')
set(haxes(2),'YColor','k')
set(haxes(1),'DataAspectRatio',[1,1,1])
set(haxes(2),'DataAspectRatio',[180,1,1])
% set(haxes(1),'square')
% set(haxes(2),'square')
[haxes,hline1,hline2]=plotyy(SortedThetas(:,2),MeanVar(:,1)-MeanVar(:,2)-5,SortedThetas(:,2),Imaginaries(:,1));
set(haxes(1),'FontSize',16)
set(haxes(2),'FontSize',16)
set(haxes(1),'XTick',0:20:180)
set(haxes(2),'XTick',0:20:180)
set(haxes(1),'YLim',[0,180])
set(haxes(1),'YTick',0:20:180)
set(haxes(2),'YLim',[0,1])
set(haxes(2),'YTick',0:.2:1)
set(hline1,'LineStyle','-')
set(hline1,'LineWidth',1.5)
set(hline1,'Color','k')
set(haxes(1),'YColor','k')
set(hline2,'LineStyle','-')
set(hline2,'LineWidth',.5)
set(hline2,'Color','k')
set(haxes(2),'YColor','k')
set(haxes(1),'DataAspectRatio',[1,1,1])
set(haxes(2),'DataAspectRatio',[180,1,1])
% set(haxes(1),'square')
% set(haxes(2),'square')
hold off;

% PRINT DATA FOR REGRESSION COMPUTATION IN EXCEL
save('rsst.txt','rSquare','-ASCII')
%     plot(SortedThetas(:,2),MeanVar(:,1)+MeanVar(:,2),'k-');
%     hold on;
%     plotyy(SortedThetas(:,2),MeanVar(:,1),SortedThetas(:,2),Imaginaries(:,2));
%     plotyy(SortedThetas(:,2),MeanVar(:,1)-MeanVar(:,2),SortedThetas(:,2),Imaginaries(:,1));
    % 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   Additional details on these methods
% %
% %
% %To draw a Random vectors from an N-dimensional multivariate normal (= Gaussian)
% %distribution w/ mean mu and covariance matrix sigma: 
% % (accomplished by the mgd.m function)
% %1- compute the cholesky decomposition (matrix square root) of sigma - see
% %for details
% %2 - generate z  = a vector whose components are N independent stanandard
% %normal variates (generate with the Box Mullter
% %transformation: http://en.wikipedia.org/wiki/Box-Muller_transform)
% %
% % A Box-Muller transform (by George Edward Pelham Box and Mervin Edgar Muller 1958)[1] is a method of generating pairs of independent standard normally distributed (zero expectation, unit variance) random numbers, given a source of uniformly distributed random numbers
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
