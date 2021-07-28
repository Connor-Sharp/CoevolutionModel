# CoevolutionModel
Individual and frequency based models for host symbiont coevolution

Coevo.m
Frequency based evolutionary model as used in Sharp and Foster

Outputs include the mean values across the population for
[Host cooperation, microbial cooperation, variance in microbial cooperation, Host Control]
Example params as used in figure 1a
coevo(2,2,0.5,0.02,0,gens,1*10^-9,0.05,.5,1,0.1,1,1*10^-6)
AB = x = Benefit to A of cooperation from B
BA = y = Benefit to B of cooperation from A
Ropt = R = Relatedness
fcost = f = indirect cost of host control
Cmax = cmax = Max host control
gen  = Host Generations
amigrate = migration of hosts into the system
bmigrate = migration of microbes into the system
stdevb = Stdard deviation of starting distribution for cooperation
stdebpa = Stadard deviation of starting distribution for control
gcost = Direct cost of control
gen2 = Within Host generations of microbes
littleB = With host generation migration of microbes into the system 

coevoIndi.m
Individual based implementation of Sharp Foster 2021.
gens = Host Generations
bmig = Migration of microbes into the system per host generation
amig = Migration of Host into the system per host generation
HostSize = Number of Hosts to model (e.g. 1*10^4)
gutSize = Number of microbes model in EACH host (e.g. 1000)
controlEffect = f = indirect cost of control
controlCost = g = direct cost of control.
