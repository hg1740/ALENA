classdef Pbar < matlab.mixin.SetGet
    %PBAR Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id       = [];
        si       = [];
        type     = [];
        mat      = [];
        a        = [];
        i        = [];
        j        = [];
        kShear   = [];
        rhoNs    = [];
        startPnt = zeros(4,2); % Stress recovery coefficients
        part     = [];
    end
    
    methods
        
        function obj = Pbar(~)
            
        end
                
        function f = genNastranDeck(obj)
            f = {};
            f = cat(1,f,sprintf('GRID*   %16.0f%16.0f%16.9E%16.9E',obj.id,obj.cs,obj.coord(1),obj.coord(2)));
            f = cat(1,f,sprintf('*       %16.9E%16.0f%16.0f%16.0f',obj.coord(3),obj.cd,obj.ps,obj.seid));
        end
        
        function f = genNeocassDeck(obj)
            id  = sprintf('%-8d', obj.id);
            mat = sprintf('%-8d', obj.mat);
            a   = sprintf('%-8s', num2_8ch(obj.a));
            I1  = sprintf('%-8s', num2_8ch(obj.i(1)));
            I2  = sprintf('%-8s', num2_8ch(obj.i(2)));
            I12 = sprintf('%-8s', num2_8ch(obj.i(3)));
            j   = sprintf('%-8s', num2_8ch(obj.j));
            
            if isempty(obj.rhoNs) || abs(obj.rhoNs)<1e-9
                rho = sprintf('%.1f', 0);
            else
                rho = sprintf('%-8g', obj.rhoNs);
            end
            
            f{1,1} = sprintf('PBAR    %s%s%s%s%s%s%s        ', id, mat, a, I1, I2, j, rho);
            
            tmp = [];
            if isempty(obj.startPnt)
                f{2} = sprintf('        0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     ');
            else
                for n = 1:size(obj.startPnt,1)
                    c1  = sprintf('%-8s', num2_8ch(obj.startPnt(n,1)));
                    c2  = sprintf('%-8s', num2_8ch(obj.startPnt(n,2)));
                    tmp = cat(2,tmp,sprintf('%s%s', c1, c2));
                end
                f{2} = sprintf('%s%s',blanks(8),tmp);
            end
            f{3} = sprintf('%s%s',blanks(24),I12);
        end
    end
end
