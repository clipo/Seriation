# Seriation #

Algorithms, implementations, manuscripts, and test cases for iterative deterministic frequency seriation done by agglomeration.

These algorithms are intended to automate and extend the initial, mainly manual implementation of IDSS in Lipo (2001).



## Directory Structure ##

* analysis
> Contains R, mathematica, and other scripts aimed at analyzing the basic problem and trial algorithms,
> before being implemented in software.  Also will contain benchmarks and performance analysis of finished
> software versions, to document scaling of the algorithms and optimizations.   
* doc
> Software documentation, installation instructions, FAQs, todo lists
* java
> Main root directory for a Java library implementing the seriation algorithms, for embedding in many types of 
> systems.  May also contain simple UI implementations (CLI, web app) to run seriations.  
* manuscripts
> Directories containing journal manuscripts describing the algorithm or applying it to data.  Each directory
> represents a separate paper and is self-contained.  
* output
> Directory for writing finished seriation solutions for test data sets.  This provides common configuration for 
> software implementations and system tests.
* perl
> Perl implementation of the seriation algorithm.  This serves as our initial testing ground and way of trying 
> ideas for optimizations, and it will be useful long-term for small or "normal" data sets.  
* testdata
> Repository of real and idealized data sets used to test the algorithm and implementations.  



### Authors ###

Carl P. Lipo, California State University at Long Beach
Mark E. Madsen, University of Washington, Seattle



(1)  Lipo, Carl P.  2001  __Science, Style and the Study of Community Structure: An Example from the Central
1317 Mississippi River Valley.__ British Archaeological Reports, International Series, no. 918.  Oxford.  
