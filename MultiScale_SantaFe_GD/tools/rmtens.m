function [SOIT,SSOIT,FOSIT,FOSPT,FODPT,FOAT1,FOAT2,SONT,SONTsqrt] = rmtens(struhyp,ntens)

    %  Tensorial operators (en notación de Voight)

    %global SOIT SSOIT FOSIT FOSPT FODPT FOAT1 FOAT2 SONT struhyp ntens

    % Second order identity tensor
    % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    SOIT = [1 1 1 0 0 0]';

    % Square second order identity tensor
    % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    SSOIT = SOIT*SOIT.';

    % Fourth order symetric indentity tensor
    % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    FOSIT = eye(6);

    % Fourth order spheric proyection tensor
    % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    FOSPT = 1/3*(SOIT*SOIT.');

    % Fourth orden deviatoric proyection tensor
    % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    FODPT = FOSIT-FOSPT;

    % Fourth orden auxiliar tensor 1
    % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    FOAT1 = [eye(3),zeros(3);zeros(3),0.5*eye(3)];

    % Fourth orden auxiliar tensor 2
    % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    FOAT2 = FOAT1*FODPT;

    % Second order norm tensor
    % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
    SONT = [eye(3),zeros(3);zeros(3),2*eye(3)];
    
    if struhyp==1||struhyp==2||struhyp == 4
        SOIT  = SOIT(1:ntens);
        SSOIT = SSOIT(1:ntens,1:ntens);
        FOSIT = FOSIT(1:ntens,1:ntens);
        FOSPT = FOSPT(1:ntens,1:ntens);
        FODPT = FODPT(1:ntens,1:ntens);
        FOAT1 = FOAT1(1:ntens,1:ntens);
        FOAT2 = FOAT2(1:ntens,1:ntens);
        SONT  = SONT (1:ntens,1:ntens);
    end

end