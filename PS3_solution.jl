include("Flux.jl");

#Parameters
enzyme = 0.01/1000; #steady-state enzyme concentration [mmol/gDW]
kcat_v1 = 203 	#[1/s]
kcat_v2 = 34.5 	#[1/s]
kcat_v3 = 249 	#[1/s]
kcat_v4 = 88.1 	#[1/s]
kcat_v5 = 13.7 	#[1/s]

Km_v1_asp = 0.15 	#[mM]
Km_v1_cit = 0.056 	#[mM]
Km_v1_ATP = 0.051 	#[mM]
Km_v2 = 3  			#[mM]
Km_v3 = 1.4 		#[mM]
Km_v4_carb_p = 0.13 #[mM]
Km_v4_orn = 0.36 	#[mM]
Km_v5 = 0.0044 		#[mM]

cit = 2.7e1 	#[mM]
asp = 1.49e1 	#[mM]
arg_s = 1 		#[mM]
fum = 4.85e-1 	#[mM]
arg = 2.18e1 	#[mM]
urea = 1 		#[mM]
orn = 4.49 		#[mM]
carb_p = 1 		#[mM]
ATP = 2.25 		#[mM]

#Enzyme velocities
v1 = (kcat_v1*enzyme)*(asp/(Km_v1_asp + asp))*(ATP/(Km_v1_ATP+ATP)); #6.3.4.5 [mmol/gDW s]
v2 = (kcat_v2*enzyme); #4.3.2.1 [mmol/gDW s]
v3 = (kcat_v3*enzyme)*(arg/(Km_v3+arg)); #3.5.3.1 [mmol/gDW s]
v4 = (kcat_v4*enzyme)*(orn/(Km_v4_orn+orn)); #2.1.3.3 [mmol/gDW s]
v5 = (kcat_v5*enzyme)*(arg/(Km_v5+arg)); #1.14.13.39 [mmol/gDW s]
v6 = (kcat_v5*enzyme); #1.14.13.39 [mmol/gDW s]
b1 = 10/3600; #Carb-P input [mmol/gDW s]
b2 = (10/3600)*((1.49*10^(-2))/(1.49*10^(-2)+1*10^(-2))); #Asp input [mmol/gDW s]
b3 = (10/3600)*((4.85*10^(-4))/(4.85*10^(-4)+5.3*10^(-3))); #Fum output [mmol/gDW s]
b4 = 10/3600; #Urea output [mmol/gDW s]
m1 = (10/3600)*((4.23*10^(-5))/(4.23*10^(-5)+6.46*10^(-5))); #AMP output [mmol/gDW s]
m2 = (10/3600)*((4.67*10^(-3))/(4.67*10^(-3)+3*10^(-5))); #ATP input [mmol/gDW s]
m3 = 10/3600; #P output [mmol/gDW s]
m4 = 10/3600; #PP output [mmol/gDW s]
m5 = 10/3600; #Water input [mmol/gDW s]

#Stochiometry
#  [v1  v2  v3  v4  v5  v6  b1  b2  b3  b4  m1  m2  m3  m4  m5]
S=[[-1   0   0   0   0   0   0   1   0   0   0   0   0   0   0];  #Asp
   [ 1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0];  #Arg-s
   [ 0   1   0   0   0   0   0   0  -1   0   0   0   0   0   0];  #Fum
   [ 0   1  -1   0  -1   1   0   0   0   0   0   0   0   0   0];  #Arg
   [ 0   0   1   0   0   0   0   0   0  -1   0   0   0   0   0];  #Urea
   [ 0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0];  #Orn
   [ 0   0   0  -1   0   0   1   0   0   0   0   0   0   0   0];  #Carb-P
   [-1   0   0   1   1  -1   0   0   0   0   0   0   0   0   0];  #Cit
   [ 1   0   0   0   0   0   0   0   0   0  -1   0   0   0   0];  #AMP
   [-1   0   0   0   0   0   0   0   0   0   0   1   0   0   0];  #ATP
   [ 0   0   0   1   0   0   0   0   0   0   0   0  -1   0   0];  #P
   [ 1   0   0   0   0   0   0   0   0   0   0   0   0  -1   0];  #PP
   [ 0   0   1   0  -1   1   0   0   0   0   0   0   0   0  -1];  #Water
   [ 0   0   0   0   1  -1   0   0   0   0   0   0   0   0   0.0]]; #Amm

#Bounds array
default_bounds_array = [[0 v1];
                        [0 v2];
                        [0 v3];
                        [0 v4];
                        [0 v5];
						[0 v6];
                        [0 b1];
                        [0 b2];
                        [0 b3];
                        [0 b4];
                        [0 m1];
                        [0 m2];
                        [0 m3];
                        [0 m4];
                        [0 m5]];

species_bounds_array = zeros(Float64, 14, 2);

objective_function_array = [0.0; #v1
                  0.0; #v2
                  0.0; #v3
                  0.0; #v4
                  0.0; #v5
                  0.0; #v6
                  0.0; #b1 
                  0.0; #b2 
                  0.0; #b3 
                 -1.0; #b4 
                  0.0; #m1 
                  0.0; #m2
                  0.0; #m3 
                  0.0; #m4 
                  0.0;] #m5 
              

#Flux calculation
flux, n1, n2, n3, n4, n5 = calculate_optimal_flux_distribution(S,default_bounds_array,species_bounds_array,objective_function_array);
flux*-3600 #Flux in [mmol/gDW h]