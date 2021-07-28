function [meanCoopA, meanCoopB, meancontrol, coopBVar] = coevoIndi(gens, bmig, amig,pmig,  HostSize, gutSize, controlEffect,controlCost)
 %coevoIndi  Individual based implementation of Sharp Foster 2021.
 %gens = Host Generations
 %bmig = Migration of microbes into the system per host generation
 %amig = Migration of Host into the system per host generation
 %HostSize = Number of Hosts to model (e.g. 1*10^4)
 %gutSize = Number of microbes model in EACH host (e.g. 1000)
 %controlEffect = f = indirect cost of control
 %controlCost = g = direct cost of control.
        generations = gens;
        controlMax = 3;
        %Set up populations
        %Make normal ditribution for coopA, controlA and coopB
        pd = makedist('normal', 'mu', 0., 'sigma',.5);
        %CoopA and coopB distribution
        t = truncate(pd,0,1);
        %controlA distribution
        z=truncate(pd, 0, controlMax);
        
        %Preallocate variables
        meancontrol = zeros(generations, 1);
        meanCoopB = zeros(generations, 1);
        meanCoopA = zeros(generations, 1);
        coopBVar = zeros(generations, 1);
        coopTemp = zeros(HostSize, 1);
        controlTemp = zeros(HostSize, 1);
        fitnessA=zeros(HostSize,1);
        fitnessB=zeros(HostSize, gutSize);
        pathB = zeros(HostSize,gutSize);
       

        
        
        %Initial pops are also used as environmental pools
        coopBTemp = random(t,HostSize, gutSize);
        oneGut = random(t,1,gutSize);
        coopBStart = repmat(oneGut,HostSize,1);
        coopAStart = random(t,HostSize,1);
        controlAStart = random(z,HostSize,1);
        coopA = coopAStart;
        coopB = coopBStart;
        controlA = controlAStart;

        %Get the stats of the starting pops
        meancontrol(1) = mean(controlA);
        meanCoopB(1)  = mean(mean(coopB,2));
        meanCoopA(1) = mean(coopA);
        coopBVar(1) = mean(var(coopB,0,2));

        

        %other model Parameters
        AB=5;
        BA=5;
        controlEffect=controlEffect;
        controlCost=controlCost;
        migA = amig;
        migB = bmig;
        migP = pmig;

        for gen = 1:generations

            %Host fitness

            gcontrol = exp(controlA .* coopB);
            changeB = sum(exp(controlA .* coopB),2);
            grelcontrol = gcontrol./changeB;
            grelcontrol2 = grelcontrol.*exp(-controlEffect*controlA);
            
            pffa = grelcontrol2.*coopA*AB;
            
            gmean = sum(pffa.*coopB,2);


            %Calculate the pathfactor (How badly a host is affected by
            %pathogens)
            immunity = controlA/controlMax;
            pathogens = sum(pathB,2)/gutSize;
            pathfactor = exp(5*pathogens.*immunity).*exp(-5*pathogens);
            
            
            fitnessA = ((1-coopA) + (gmean*BA) - controlCost*controlA).*pathfactor;
            fitnessA(fitnessA < 0) = 0;
            fitnessA = fitnessA./sum(fitnessA);
            
            %make new populations based on fitness values
            weights = randsample(linspace(1,HostSize,HostSize), HostSize, true,  fitnessA);
            controlTemp = controlA(weights);
            coopTemp = coopA(weights);
            coopBTemp = coopB(weights,:);

            coopB = coopBTemp;
            controlA = controlTemp;
            coopA = coopTemp;

            %microbesl fitness
            changeB = sum(exp(controlA .* coopB)/gutSize,2);
            expCJ = exp(-controlEffect*controlA);
            gcontrolAB = exp(controlA .* coopB)/gutSize;
            grelcontrolAB = gcontrolAB./changeB;
            grelcontrol2AB = grelcontrolAB .* expCJ;
            pffb = sum((coopB.*grelcontrol2AB),2)*BA;
            pffbAB = pffb * AB;
            gcontrolB = exp(controlA .* coopB);

            gcontrolBrel = gcontrolB./changeB;

            gcontrolBrel2 = gcontrolBrel .* expCJ;

            gcontrolcoop = gcontrolBrel2.*coopA;
            fitnessB = (1 - coopB) + (gcontrolcoop.*pffbAB);
            fitnessB = fitnessB./sum(fitnessB,2);

            for i = 1:HostSize
                index = randsample(gutSize, gutSize, true,  fitnessB(i,:));
                coopB(i,:) = coopB(i,index);
                pathB(i,:) = pathB(i,index);
            end

            meancontrol(gen+1) = mean(controlA);
            meanCoopB(gen+1)  = mean(mean(coopB,2));
            meanCoopA(gen+1) = mean(coopA);
            coopBVar(gen+1) = mean(var(coopB,0,2));

            %add migration in

            % Host migration
            changeME = randi([1,HostSize],HostSize*migA,1);
            coopA(changeME) = coopAStart(changeME);
            controlA(changeME) = controlAStart(changeME);
            coopA(changeME) = 0.1*coopAStart(changeME);
            controlA(changeME) = 0.1*controlAStart(changeME);

            % microbesl migration
            changeHost = randi([1,HostSize],gutSize*migB,1);
            changegut = randi([1,gutSize],gutSize*migB,1);
            coopB(changeHost,changegut) = coopBStart(changeHost, changegut);
            
            %Pathogen migration - almost identical to normal migration but
            %uses pmig and updates the pathogen matrix aswell
            changeSize = round(gutSize*(migP*migB));
            
            %make new populations based on fitness values
            if changeSize >= 1

                changeHost = randi([1,HostSize],changeSize,1);
                changegut = randi([1,gutSize],changeSize,1);
                coopB(changeHost,changegut) = 0;
                pathB(changeHost, changegut) = 1;
            end
  
            %Mix microbes between host in the generation.
            %This is similar to mixing microbes in a community.
            coopBTemp = coopB;
            coopB = reshape(coopBTemp(randperm(gutSize*HostSize)), HostSize, gutSize );

        end

end








