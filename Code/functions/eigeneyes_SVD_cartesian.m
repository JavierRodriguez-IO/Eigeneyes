function [RMS,kp5,coefs_correlation,polynomials_third_order,Reconstructed_eyes,signals,tv,mn,PC,cum_explained] = eigeneyes_SVD_cartesian(data,x_sampling,y_sampling,var_max,Num_lenses,fixed_number_of_eigen, age_vector)
%
% function [RMS,kp5,coefs_correlation,polynomials_third_order,Reconstructed_eyes,signals,tv,mn,PC] = eigeneyes_SVD_cartesian(data,x_sampling,y_sampling,var_max,Num_lenses,fixed_number_of_eigen, age_vector)
%
% example call: % TEST CODE
%               % [RMS,kp5,coefs_correlation,polynomials_third_order,Reconstructed_eyes,signals,tv,mn,PC]= eigeneyes_SVD_cartesian(dataZ_eigeneyes_pre,dataX_eigeneyes_pre,dataY_eigeneyes_pre,var_max,Num_eyes,age_vector);
%
% INPUT:
% data:                  distance of [ZAC;ZPC;ZAL;ZPL;ZRE].
% x_sampling:            X-axis sampling.
% y_sampling:            Y-axis sampling.
% var_max:               max variance in the eigeneye representation.
% Num_lenses:            number of eyes of the database.
% fixed_number_of_eigen: number of eigeneyes used in the reconstruction.
% age_vector:            age vector of the eyes in the database.
%
% OUTPUT:
% RMS of the lenses reconstructed with kp5 eigenvectors;
% kp5: Number of eigenvectors used in the reconstruction
% coefs_correlation: CORRELATION COEFFICIENTS BETWEEN THE FIRST kp5 EIGENLENSES AND AGE USING LINEAR MODELS
% polynomials_third_order: MATRIX WHERE EVERY ROW IS A EIGEN COEF, AND
% EVERY COLUMN THE COEFFICIENT OF THE POLYNOMIAL MODEL THAT REPRSENT COEFFICIENT F(AGE)
% Reconstructed_eyes: RECONSTRUCTION OF THE LENSES WITH KP5 EIGENLENSES AND RMS CALCULATION
% signals: Coefs in the new basis (eigeneyes) for every lens
% tv: percentaje of variance explained with the kp5 coefficients (sum of
% eigenvalues)
% mn: mean lens
% PC: Principal components
%
% This script applies the SVD and calculate the eigeneyes and the
% cumulative explained variance with these eigeneyes.

% Eigeneyes: A statistical model of whole-eye geometry and its applications
% Copyright (C) 2026 Estela Vadillo, Álvaro de la Peña,
% Javier Rodríguez-Sánchez, Alberto de Castro, Susana Marcos,
% Eduardo Martínez-Enríquez
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Variables definition to
RMS                = zeros(1,Num_lenses);
Reconstructed_eyes = zeros(length(x_sampling),Num_lenses);
% kp5 =6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBTAIN THE EIGENLENSES (USING SVD OR PCA). SVD IS FASTER,  PCA IS MORE
% INTUITIVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
% SVD
%%%%%%%%%%%%%
% mean lens
data_original = data;
mn = mean(data,2);
% data= original lenses data- mean lens (residual for each lens)
data = data-repmat(mn,1,Num_lenses);
Y = data'/sqrt(Num_lenses-1);
[u,S,PC] = svd(Y);
S = diag(S);
eigval = S.*S;
% proyection of the residual data (data) in the eigenlenses (PC)--->
% signals (coefs in the new basis)--> rotation of the data
signals = PC'*data;

% % % % % %%%%%%%%%%%%%
% % % % % % PCA
% % % % % %%%%%%%%%%%%%
% % % % %
% % % % % covariance=1/(Num_lenses-1)*data*data';
% % % % % % obtain eigenvalue & eigenvector
% % % % % [PC,V]=eig(covariance);
% % % % % % % % % % % sort the variances in decreasing order
% % % % % V=real(V);
% % % % % eigval=diag(V);
% % % % % [junk,rindices]=sort(-1*eigval);
% % % % % eigval=eigval(rindices);
% % % % % PC=PC(:,rindices);
% % % % %
% % % % %  % PROYECTION
% % % % % signals=PC'*data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECT THE NUMBER OF EIGENLENSES USED IN THE RECONSTRUCTION AND INITALIZE
% VARIABLES WITH THE CORRESPONIDING SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

explained = eigval / sum(eigval);
cum_explained = cumsum(explained);
eigsum = sum(eigval);
csum = 0;
if  fixed_number_of_eigen==0
    for i=1:length(data)
        csum = csum+eigval(i);
        tv   = csum/eigsum;
        if tv>var_max
            kp5 = i;
            break
        end
    end
else
    kp5 = fixed_number_of_eigen;
    for i=1:kp5
        csum = csum+eigval(i);
        tv   = csum/eigsum;
    end
end

% kp5= 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBTAIN CORRELATION COEFFICIENTS BETWEEN THE FIRST kp5 EIGENLENSES AND AGE USING LINEAR MODELS.
% ALSO OBTAIN EACH OF THE FIRST KP5 COEFFICIENTS F(AGE) USING NON-LINEAR MODELS (E.G., POLYNOMIAL MODELS OF ORDER 3)
% (USEFUL, FOR INSTANCE, TO RECONSTRUCT LENSES OF A GIVEN AGE). MATRIX "polynomials_third_order"
% WHERE EVERY ROW IS A EIGEN COEF, AND EVERY COLUMN THE COEFFICIENT OF THE POLYNOMIAL MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coefs_correlation = zeros(kp5,1);
polynomials_third_order = zeros(kp5,4);
% for j=1:kp5
%     coefs = signals(j,:);
%     [r,p] = corrcoef(age_vector,coefs);
%     coefs_correlation(j) = r(1,2);
%     % IN SIGNALS WE HAVE THE COEFFICIENTS OF THE RESIDUAL DATA IN THE
%     % EIGENBASIS (SEE ABOVE)
%     P = polyfit(age_vector', signals(j,:),3);
%     polynomials_third_order(j,:) = P;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECONSTRUCTION WITH KP5 EIGENLENSES AND RMS CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Num_lenses
    lens_aux = mn;
    for j=1:kp5
        %    lens reconstruction=mean_lens+new_eigen_basis (rotation)*eigen_coeficients
        lens_aux = mn+PC(:,[1:j])*signals((1:j),i);
    end
    % RMS between original lens (residual_data+mean lens) and reconstructed
    % from kp5 eigenlenses
    RMS(i) = rms(data(:,i)+mn-lens_aux);
    Reconstructed_eyes(:,i) = lens_aux;
end

end