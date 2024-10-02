function p = compare_correlation_coefficients(r1,r2,n1,n2)

% This function compare if two correlation coefficients are significantly
% different.
% The correlation coefficients were tansfered to z scores using fisher's r
% to z transformation. 
% ref: http://core.ecu.edu/psyc/wuenschk/docs30/CompareCorrCoeff.pdf
%--------------------------------------------------------------------------
% Inputs: (1) r1: correlation coefficient of the first correlation (2) r2:
% correlation coefficient of the second correlation (3) n1: number of
% samples used to compute the first correlation (4) n2: number of samples
% used to compute the second correlation
%--------------------------------------------------------------------------
% Output: (1) p: p value, the probability that H0 (the correlation
% coefficiets are not different) is correct
%--------------------------------------------------------------------------
% Example :
% x = rand(20,1); 
% y1= x+rand(20,1)*0.05;
% y2= x+rand(20,1)*0.5;
% r1=corr(x,y1);
% r1=corr(x,y2);
% p = compare_correlation_coefficients(r1,r2,length(x),length(x));
%--------------------------------------------------------------------------
% Copyright (c) 2013, Sisi Ma
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%--------------------------------------------------------------------------

t_r1 = 0.5*log((1+r1)/(1-r1));
t_r2 = 0.5*log((1+r2)/(1-r2));
z = (t_r1-t_r2)/sqrt(1/(n1-3)+1/(n2-3));
p = (1-normcdf(abs(z),0,1))*2;
end