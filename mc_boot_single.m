function [p_value] = mc_boot_single(data,tick)

open=data(:,4);
if(min(min(open))<=0) % check for negative/zero values due to continuous backadjusted prices
open=(abs(min(min(open)))+tick).+open; % if they exist, compensate
endif

%***********************************************************************************************
% Create the position vector for the test in question - code to be changed according to the test
%***********************************************************************************************

long = ( data(:,142) > 0 );
short = ( data(:,142) < 0 ).*-1; short( short == -0 ) = 0;
pos_vector = long .+ short;

%***********************************************************************************************
% end of code for the position vector 
%***********************************************************************************************

log_return_vector = log10 ( open ./ shift(open , 1 ) ) ; log_return_vector(1,1) = 0 ;

% shift the position vector to "agree" with above return vector
pos_vector = shift(pos_vector,2); pos_vector(1:2,1)=0; 

% truncate log_return_vector & pos_vector to account for above shift and any lead-in time to allow indicator(s) to settle
final_pos_vector = pos_vector(50:end,1);
detrended_log_return_vector = center(log_return_vector(50:end,1));

p_value = mc_MTboot_single(detrended_log_return_vector,final_pos_vector,5000);

endfunction



 



