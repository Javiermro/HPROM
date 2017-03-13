function [inclambda] = factorlambda (Fext,du);

inclambda = (Fext.^2 + du.^2) + inclambda;
