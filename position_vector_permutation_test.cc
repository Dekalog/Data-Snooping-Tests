#include <octave/oct.h>
#include <octave/dRowVector.h>
#include <octave/dMatrix.h>
#include <math.h>
#include <algorithm>
#include "MersenneTwister.h"
     
DEFUN_DLD ( position_vector_permutation_test, args, nargout, 
"-*- texinfo -*-\n\
@deftypefn {Function File} {} position_vector_permutation_test (@var{returns,position_vector_matrix,iters})\n\
This function takes as inputs a vector of returns, a position vector matrix the columns of which are\n\
the individual position vectors that are to be part of the test and the number of iterations for the test.\n\
If the matrix consists only of one position vector the test performed is the single position vector permuation\n\
test. If there are two or more position vectors the test performed is the data mining bias adjusted permutation\n\
test. The null hypothesis that is tested is that the position vector(s) give results that are no better than those\n\
a random ordering of equivalent net positions would give. The alternative hypothesis is that the position vector(s) give\n\
results that are better than a random ordering would give. The function returns a p-value(s) for the position vector(s)\n\
being tested. If small enough, we can reject the null hypothesis and accept the alternative hypothesis. The returned\n\
p-values are expressed as percentages, i.e 0.05 means a 5% significance value. The test statistic is the sum of returns.\n\
*IMPORTANT NOTE* Apart from basic input argument checks to avoid error messages and seg faults this function does not\n\
check that the returns vector and position vector matrix are properly alligned etc. with regard to accuracy of test results\n\
and makes no assumptions about this. It is the user's responsibility to ensure that the position vector(s) are offset to\n\
match the returns for the particular test in question.\n\
@end deftypefn" )

{
octave_value_list retval_list ; // the return value(s) are p value(s)
int nargin = args.length () ;

// check the input arguments
if ( nargin != 3 )
   {
   error ("Insufficient arguments. Inputs are a returns vector, a position vector matrix and a scalar value for iters.") ;
   return retval_list ;
   }

if ( args(0).length () != args(1).rows () )
   {
   error ("Dimensions mismatch. Length of returns vector should == No.rows of position_vector_matrix.") ;
   return retval_list ;
   }
   
if ( args(2).rows() > 1 )
   {
   error ("Invalid 3rd argument. This should be a scalar value for the number of permutations to perform.") ;
   return retval_list ;
   }
   
if ( args(2).length () != args(2).rows () )
   {
   error ("Invalid 3rd argument. This should be a scalar value for the number of permutations to perform.") ;
   return retval_list ;
   }   

if ( error_state )
   {
   error ("Error state. See usage") ;
   return retval_list ;
   }
// end of input checking

//disable_warning ( "Octave:broadcast" ) ;
set_warning_state ( "Octave:broadcast" , "off" ) ;

RowVector returns = args(0).row_vector_value () ; // Returns vector input
Matrix position_vector = args(1).matrix_value () ; // Positions vector input

// compute the return(s) of the candidate model(s) 
Matrix model_returns ( 1 , position_vector.cols() ) ;
model_returns = returns * position_vector ; // matrix multiplication

// normalised returns multiplier
Matrix normalised_multiplier ( 1 , position_vector.cols() ) ;
normalised_multiplier.fill ( 0.0 ) ;
   
if ( position_vector.cols() > 1 ) // more than one position vector
   {
   for ( octave_idx_type ii (0) ; ii < position_vector.cols() ; ii++ )
       { 
	 for ( octave_idx_type jj (0) ; jj < position_vector.rows() ; jj++ )
	     {
	     normalised_multiplier( 0 , ii ) += fabs( position_vector( jj , ii ) ) ;  
	     } // end of jj loop     
	 normalised_multiplier( 0 , ii ) = sqrt( position_vector.rows() ) / sqrt( normalised_multiplier( 0 , ii ) ) ; 
       } // end of ii loop   
   } // end of if ( position_vector.cols() > 1 )   
   else // only one position vector
   {
   normalised_multiplier( 0 , 0 ) = 1.0 ;
   }

// use this normalised_multiplier to adjust the model returns
model_returns = product( model_returns , normalised_multiplier ) ;  

// Now do the Monte-Carlo replications
int temp ; 
int k1 ;
int k2 ;
Matrix trial_returns ( 1 , position_vector.cols() ) ;
double max_trial_return ;
Matrix counts ( 1 , position_vector.cols() ) ;
counts.fill ( 0.0 ) ;
MTRand mtrand1 ; // Declare the Mersenne Twister Class - will seed from system time

 for ( octave_idx_type ii (0) ; ii < args(2).int_value() ; ii++ ) // args(2).int_value() is the no. of permutaions to perform 
     {
       trial_returns.fill( 0.0 ) ; // set trial returns to zero

       k1 = args(0).length () - 1 ; // initialise prior to shuffling the returns vector

       while ( k1 > 0 ) // While at least 2 left to shuffle
             {          
             k2 = mtrand1.randInt ( k1 ) ; // Pick an int from 0 through k1 

               if (k2 > k1)     // check if random vector index no. k2 is > than max vector index - should never happen 
                  {
                  k2 = k1 - 1 ; // But this is cheap insurance against disaster if it does happen
                  }

             temp = returns ( k1 ) ;           // allocate the last vector index content to temp
             returns ( k1 ) = returns ( k2 ) ; // allocate random pick to the last vector index
             returns ( k2 ) = temp ;           // allocate temp content to old vector index position of random pick
             k1 = k1 - 1 ;                     // count down 
             }                                 // Shuffling is complete when this while loop exits
             
      // Compute returns for this random shuffle
      trial_returns = returns * position_vector ;
      
      // and now normalise these returns
      trial_returns = product( trial_returns , normalised_multiplier ) ;
      
      max_trial_return = *std::max_element( &trial_returns( 0,0 ), &trial_returns( 0,trial_returns.cols() ) ) ;

        for ( octave_idx_type jj (0) ; jj < trial_returns.cols() ; jj++ )
	    {
            if ( max_trial_return >= model_returns(0,jj) ) // If the best random system(s) beats the candidate system(s)
               {
               counts(0,jj) += 1 ; // Count it
               }
	    }
	    
     } // end of main ii loop

// now write out the results     
Matrix p_values ( 1 , position_vector.cols() ) ;
Matrix iters = args(2).matrix_value () ;   
p_values = quotient( counts , iters ) ; // count divided by number of permutations
 
retval_list(0) = p_values ;
 
set_warning_state ( "Octave:broadcast" , "on" ) ;

return retval_list ; // Return the output to Octave

} // end of function