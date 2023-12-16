###### DATA_AVAILABILITY ######

    This repository contains 1 file:
    
 
    (1) One ".cpp" file for the estimate of the rescue/extinction probability:
	
			C++ source code "rescue-L.cpp"
	
			Requires the GSL scientific library:
		
			Ubuntu systems installation through terminal:
			
	        (C++)	~ sudo apt-get install build-essential
			(GSL)	~ sudo apt-get install gsl-bin
			(GSL)	~ sudo apt-get install libgsl-dev
			
			Ubuntu systems compilation through terminal:
		
			~ c++ -O3 rescue-L.cpp -o [executable_name] -lm -lgsl -lgslcblas
	
			After compilation, we will get an executable.

                       To run the code, we must provide it with input data.
                        
                       The input data are:
                          
                       i- sequence size
                       ii- carrying capacity
                       iii - mutation probability
		       iv - number of independent landscapes
		       v - maximum time (usually before this maximum as a consequence of extinction, or because the population reached the global optimum. The accessibility of the global optimum is extremely high. Because in the genotypic FGM there can exist other local optimum, in some few simulations the population can use this maximum time. So, usually, we make this time very long, around 10,000.)
                       vi - number of traits
		       vii - selection strength, which me make equal to 1
                       viii - mean value of phenotypic effects due to mutations
                       ix - initial fitness drop (the initial fitness is equal to (1 - drop_fitness)
                       x - Maximum fitness W_max

                       Here is an example of how to run the code in an Ubuntu terminal

                      ./[executable_name] 12  10000  0.005    50000    10000    5     1.0        0.4        0.37      1.5 

                      The output of the code is created in an .DAT file. The name of the file is RESCUE and also includes the values of parameter used.

                      The output contains 2 columns: value of the initial drop in fitness (parameter), the probability of extinction

                      The output was made simpler than in our version, which includes several other measurements, some not used.



    (2) One ".cpp" file for the estimate of the degree and parallelism/mean path divergence:
	
			C++ source code "parallel_evolution.cpp"
	
			Requires the GSL scientific library:
		
			Ubuntu systems installation through terminal:
			
		        (C++)	~ sudo apt-get install build-essential
			(GSL)	~ sudo apt-get install gsl-bin
			(GSL)	~ sudo apt-get install libgsl-dev
			
			Ubuntu systems compilation through terminal:
		
			~ c++ -O3 parallel_evolution.cpp -o [executable_name] -lm -lgsl -lgslcblas
	
			After compilation, we will get an executable.

                       To run the code, we must provide it with input data.
                        
                       The input data are:
                          
                       i- sequence size
                       ii- carrying capacity
                       iii - mutation probability
		       iv - number of independent landscapes
		       v - maximum time (usually before this maximum as a consequence of extinction, or because the population reached the global optimum. The accessibility of the global optimum is extremely high. Because in the genotypic FGM there can exist other local optimum, in some few simulations the population can use this maximum time. So, usually, we make this time very long, around 10,000.)
                       vi - number of traits
		       vii - selection strength, which me make equal to 1
                       viii - mean value of phenotypic effects due to mutations
                       ix - initial fitness drop (the initial fitness is equal to (1 - drop_fitness)
                       x - Maximum fitness W_max

                       Here is an example of how to run the code in an Ubuntu terminal

                      ./[executable_name] 12  10000  0.005    50000    10000    5     1.0        0.4        0.37      1.5 

                      The output of the code is created in two .DAT file. The names of the files are Parallel_MAX and Descent, and also includes the values of parameter used.

                      The file Parallel_MAX presents the values of the metrics of evolutionary pathways built storing the information about the highest fitness

                      The file Descent presents the values of the metrics of evolutionary pathways built using the line of descent

                      The output of the first file contains 8 columns: value of the initial drop in fitness (parameter), degree of parallelism, its corresponding standard error, parallelism of ending points, its corresponding standard error, mean path divergence, its corresponding standard error, accessibility of the global optimum

                      The output of the second file contains 5 columns: value of the initial drop in fitness (parameter), degree of parallelism, its corresponding standard errormean path divergence, its corresponding standard error

