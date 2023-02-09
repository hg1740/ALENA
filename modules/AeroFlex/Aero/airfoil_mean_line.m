%% airfoil_mean_line
% The following script takes a string [foil] and works out the camber angle
% 
% 
%   Author: Dario Calderon 

function [EP, Yu, Yl, xa, angle] = airfoil_mean_line(foil,plotfoil)

if isempty(str2num((cell2mat(foil))))==0
    
    TYPE = 1;       % Naca xxxx profile, case 1
    
elseif isempty(str2num((cell2mat(foil))))
    
    TYPE = 2;       % Airfoil from file, case 2
    
else
    
    error('Foil error, flatplate assumed')
    
end

STEP = 30;

switch TYPE
    
    case 1
        % NACA 4 digits
        foil = str2num(cell2mat(foil));
        if foil == 0
            foil = 0012;
        end
        m = fix(foil / 1000);	% gives first NACA-4 number
        lemma = foil - m*1000;
        p = fix(lemma/100);	  % gives second NACA-4 number
        p = p/10;
        m = m/100;
        xa = (0:1/STEP:1);
        for i=1:STEP+1
            if xa(i) < p
                a(i)=(2*m/(p^2)*(p-xa(i)));
            else
                a(i)=2*m/((1-p)^2)*(p-xa(i));
            end
        end
        angle = atan(a);
        EP = xa';
        c = 1;
        t = (foil - m*1000 - p*100)/100;
        % ONLY SETUP FOR SYMMETRIC AIRFOIL
        Yu = 5*t*c*(0.2969*sqrt(xa./c) + (-0.1260)*(xa./c) + (-0.3516)*(xa./c).^2 + 0.2843*(xa./c).^3 + (-0.1015)*(xa./c).^4);
        Yl = -5*t*c*(0.2969*sqrt(xa./c) + (-0.1260)*(xa./c) + (-0.3516)*(xa./c).^2 + 0.2843*(xa./c).^3 + (-0.1015)*(xa./c).^4);
        
    case 2
        
        
        A = load(char(foil));
        % Take the number of data points in the data file
        Nu = A(1,1); % for the upper surface
        Nl = A(1,2); % for the lower surface
        %
        % check if format is ok
        if (Nu + Nl ~= size(A,1)-1)
            errmsg = ['Airfoil file ', char(foil), ' has no upper and lower points declaration at first line or wrong values given.'];
            error(errmsg);
        end
        %
        xup = A(2:Nu+1,1);
        yup = A(2:Nu+1,2);
        xdw = A(Nu+2:end,1);
        ydw = A(Nu+2:end,2);
        [xup, index] = sort(xup);
        yup = yup(index);
        [xdw, index] = sort(xdw);
        ydw = ydw(index);
        if (Nu ~= Nl)
            % determine missing points
            X1 = setdiff(xdw, xup);
            X2 = setdiff(xup, xdw);
            yup = [yup; interp1(xup, yup, X1, 'pchip')];
            ydw = [ydw; interp1(xdw, ydw, X2, 'pchip')];
            xup = [xup; X1];
            xdw = [xdw; X2];
            [xup, index] = sort(xup);
            yup = yup(index);
            [xdw, index] = sort(xdw);
            ydw = ydw(index);
        end
        %Upper surface
        Xu = xup/(xup(end) - xup(1));
        Yu = yup/(xup(end) - xup(1));
        % Lower surface
        Xl = xdw/(xdw(end) - xdw(1));
        Yl = ydw/(xdw(end) - xdw(1));
        %
        EP = [0:1/STEP:1];
        Yu = interp1(Xu, Yu, EP, 'pchip')';
        Yl = interp1(Xl, Yl, EP, 'pchip')';
        % Mean line
        ml = 0.5.*(Yu+Yl);
        %
        %     solve least square problem
        A = [EP.^3; EP.^2; EP; ones(size(EP))]';
        sol =(A'*A)\(A'*ml);
        angle_me = atan(sol(3) + sol(2).*(2.*EP) + sol(1).*(3.*EP.^2) );
        xa = EP;
        
        for ll = 1:length(ml)
            if ll ==1
                angle(ll) = atan((ml(2)-ml(1))/(EP(2)-EP(1)));
            elseif ll ==length(ml)
                angle(ll) = atan((ml(ll)-ml(ll-1))/(EP(ll)-EP(ll-1)));
            else
                angle(ll) = 0.5*((atan((ml(ll)-ml(ll-1))/(EP(ll)-EP(ll-1))))+(atan((ml(ll+1)-ml(ll))/(EP(ll+1)-EP(ll)))));
            end
        end
        
        if plotfoil == 1
            figure(100);
            plot(EP,Yu,'.');
            hold on
            plot(EP,Yl,'.');
            plot(EP,ml,'-r.');
            plot(EP,angle,'-c.');
            plot(EP,angle_me,'-k.');
            axis equal
            hold off;
            %legend('Lower Co-ord','Lower Co-ord.','mean camber line','angle','angle tornado');
            %title(['Mean camber line and induced angle for the ' cellstr(foil) ' foil']);
        end
        
        EP = EP';
        Yu = Yu';
        Yl = Yl';
end
%
end