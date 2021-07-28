function [xamean, xbmean, varxb, pamean] = coevo(AB,BA,Ropt,fcost,Cmax,gen,amigrate,bmigrate,stdevb,stdevpa,gcost,gen2, littleB);
% coevo  Frequency based evolutionary model as used in Sharp and Foster
% (2021).
% Outputs include the mean values across the population for
% [Host cooperation, microbial cooperation, variance in microbial
% cooperation, Host Control]
%Example params as used in figure 1a
%coevo(2,2,0.5,0.02,0,gens,1*10^-9,0.05,.5,1,0.1,1,1*10^-6)
%AB = x = Benefit to A of cooperation from B
%BA = y = Benefit to B of cooperation from A
%Ropt = R = Relatedness
%fcost = f = indirect cost of host control
%Cmax = cmax = Max host control
%gen  = Host Generations
%amigrate = migration of hosts into the system
%bmigrate = migration of microbes into the system
%stdevb = Stdard deviation of starting distribution for cooperation
%stdebpa = Stadard deviation of starting distribution for control
%gcost = Direct cost of control
%gen2 = Within Host generations of microbes
%littleB = With host generation migration of microbes into the system

% Generates row vector with 11 points between 0 and 1 or cmax for coop and
% control
xa=linspace(0,1,11)';

pa=linspace(0,Cmax,11);
xb=linspace(0,1,11)';

f=normpdf(xa,0,stdevb)*normpdf(pa,0,stdevpa); 
% Normpdf Generates normal probability distributions for species A coop (x) and species A control (y) mean is 0.8 and sd is .1
% reported values are from the range defined by the vector - for x it will report the values from the
%distribution that lies between 0 and 1 on the x axis - it is finite distribution with the ends cut off. 
% f is the product of a row and a column vector which produces a 11 x 11 matrix.
fsource=f./sum(sum(f));
f=fsource;
% sum(sum(f)) returns the sum of all elements in the matrix (sum(f) is the sum of the columns)
% This makes the matrix to sum to 1 and turns things into proportions.
% f now defines the frequencies of species A with coop at x and control of y 

g=normpdf(xb,0,stdevb); gsource=g./sum(g); 
g=gsource;

%Intialise vectors for outputting results
xamean(1)= sum(sum(f')'.*xa);
pamean(1)= sum(sum(f).*pa);
xbmean(1)= sum(xb.*g);
time(1) = 1;
varxa(1)=var(sum(f'),1);
%varpa(t)=var(sum(f),pmax);
varpa(1)=var(sum(f),1);
varxb(1)=var(g',1);

gHolder = zeros(length(pa),gen2,length(pa));
gFinal = zeros(length(pa),length(pa));

% NEED TO DEFINE THE FITNESS OF EACH POSITION IN THE MATRIX IN TERMS OF THE LEVEL OF COOP AND control IN EACH SPECIES....
for t=1:gen
   for h=1:gen2
       if h==gen2
           R=Ropt;
       else
           R=0;
       end
       SumfitB = zeros(length(xb),1)';
        for k=1:length(xb)
      % For each level of coop by B, the following calculates their fitness for all levels of coop and control in A
      % At the end this is multiplied by the freq of coop and control in A to give the effects of this on B
          fitB = zeros(length(xa), length(pa));
          for j=1:length(pa)
              
              %%for each xb level
              gcontrolB = exp(pa(j).*xb(k));
              changeB = (R.*gcontrolB)+(1-R).*sum((exp(pa(j).*xb)).*g);
              gcontrolBrel = gcontrolB/changeB;
              gcontrolBrel2 = gcontrolBrel.*exp(-fcost*pa(j));
              
              %%for average xb
              gcontrolAB = (exp(pa(j).*(1-R).*xb)).*g;
              grelcontrolAB = gcontrolAB/sum(gcontrolAB);
              grelcontrol2AB = grelcontrolAB.*exp(-fcost*pa(j));
              pffb=(xb(k)*R+ (1-R).* sum(xb.*grelcontrol2AB) )*BA;
              gHolder(j,h,:) = grelcontrol2AB;
             for  i=1:length(xa)


                gcontrolcoop = gcontrolBrel2.*xa(i);

                 if h==gen2
                    fitB(i,j)=(1-xb(k))+gcontrolcoop.*pffb.*AB;
                 else
                    fitB(i,j)=(1-xb(k))+(gcontrolBrel2*AB);


                 end
             end
          end
          %Creates a matrix with the fitness effects of all levels of coop and control by species A
          %This makes the effects dependent on the actual frequency of each of the combinations of coop and control in B
          fitB = fitB.*f;
          SumfitB(k)=sum(sum(fitB));

        end
        
    gnew=g.*SumfitB';
    gnew=gnew./sum(gnew);
    g=((1-(littleB)).*gnew)+(littleB).*gsource;
    
   end


%Calculate the cooperation recieved  the host by averaging the cooperation of
%microbes at each microbial generation
for j=1:length(pa)
    for h=1:gen2
        holdingMatrix(h,:) = gHolder(j,h,:); 
    end
    gFinal(j,:) = mean(holdingMatrix,1);

end

          

   for i=1:length(xa)
 
      
      for j=1:length(pa)

         grelcontrol2 = gFinal(j,:)';

         pffa = grelcontrol2.*xa(i)*AB; %AB = y?

         gmean = sum(pffa.*xb);

       


            
         % Gives the mean level of cooperation you get back accounting for both control and your cooperativity
         if Cmax ==0
            fitAA(i,j)=((1-xa(i))+(gmean.*BA));
         else
            fitAA(i,j)=((1-xa(i))+(gmean.*BA))-(gcost.*(pa(j)/Cmax));
         end
         fnew(i,j)=f(i,j).*fitAA(i,j);
            
   
            
      end      
   end
  %Normalise
    fnew=fnew./sum(sum(fnew));
    %Add a small migration of A.
    f=((1-amigrate).*fnew)+amigrate.*fsource;
    f=f./sum(sum(f));

   
   

   
   %%%%Update all stats
   xbmean(t+1)= sum(xb.*g);
   xamean(t+1)= sum(sum(f')'.*xa);
   pamean(t+1)= sum(sum(f).*pa);
   varxb(t+1)=var(g',1);
   %%% Then add migration into microbes
   g=((1-(bmigrate)).*gnew)+(bmigrate).*gsource;


end

