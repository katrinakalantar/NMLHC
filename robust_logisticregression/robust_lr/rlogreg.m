function varargout = rlogreg(varargin)
% 
% RLOGREG - robust logistic regression with Bayesian regularisation of a Laplace prior 
%
%    RLOGREG implements a modified form of the sparse logistic regression
%    algorithm of Shevade and Keerthi [1], such that the regularisation
%    parameter is integrated out analytically, as described in [2].
%
%    References
%
%    [1] S. K. Shevade and S. S. Keerthi, "A simple and efficient algorithm
%        for gene selection using sparse logistic regression", Bioinformatics,
%        vol. 19, no. 17, pp. 2246-2253, 2003.   
%
%    [2] G. C. Cawley and N. L. C. Talbot, "Gene selection in cancer
%        classification using sparse logistic regression with Bayesian
%        regularisation", Bioinformatics (submitted), 2006.
%

%
% File        : rlogreg.m
%
% Date        : Saturday 24th March 2006
%
% Author      : Jakramate Bootkrajang
%
% Description : User documentation for MEX implementation of sparse logistic
%               regression with Bayesian regularisation using a Laplace prior.
%               Also  provides a simple veneer for tranparently compiling the
%               MEX file the first time they are used.
%
% History     : 02/10/2004 - v1.00
%               25/03/2006 - v1.10 Minor improvements to help comments etc.
%
% Copyright   : (c) Jakramate Bootkrajang 2011.
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%

% store the current working directory

cwd = pwd;

% find out name of the currently running m-file (i.e. this one!)

name = mfilename;

% find out what directory it is defined in

dir             = which(name);
dir(dir == '\') = '/';
dir             = dir(1:max(find(dir == '/')-1));

% try changing to that directory

try

   cd(dir);

catch

   % this should never happen, but just in case!

   cd(cwd);

   error(['unable to locate directory containing ''' name '.m''']);

end

% try recompiling the MEX file

try

   mex([name '.c'], '-lm');

catch

   % this may well happen happen, get back to current working directory!

   cd(cwd);

   error('unable to compile MEX version of ''%s''%s\n%s%s', name, ...
         ', please make sure your', 'MEX compiler is set up correctly', ...
         ' (try ''mex -setup'').');


end

% change back to the current working directory

cd(cwd);

% refresh the function and file system caches

rehash;

% try to invoke MEX version using the same input and output arguments

[varargout{1:nargout}] = feval(name, varargin{:});

% bye bye...

