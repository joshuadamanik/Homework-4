classdef quasi_newton_class < handle    
    properties
        H
        X_last
        g_last
    end
    
    methods
        function obj = quasi_newton_class(H0)
            obj.H = H0;
        end
        
        function p = rankone(obj, X, eps, f)
            g = grad_central_diff(X, eps, f);
            if ~isempty(obj.X_last)
                del_x = X - obj.X_last;
                del_g = g - obj.g_last;
                
                num = del_x - obj.H * del_g;
                den = del_g'*(del_x-obj.H*del_g);
                
                obj.H = obj.H + (num * num'./ den);
            end
            
            p = -obj.H*g;
            obj.X_last = X;
            obj.g_last = g;
        end
        
        function p = dfp(obj, X, eps, f)
            g = grad_central_diff(X, eps, f);
            if ~isempty(obj.X_last)
                del_x = X - obj.X_last;
                del_g = g - obj.g_last;
                
                num1 = del_x*del_x';
                den1 = del_x'*del_g;
                
                num2 = obj.H*del_g;
                den2 = del_g'*obj.H*del_g;
                
                obj.H = obj.H + num1/den1 - num2*num2'/den2;
            end
            
            p = -obj.H*g;
            obj.X_last = X;
            obj.g_last = g;
        end
        
        function p = bgfs(obj, X, eps, f)
            g = grad_central_diff(X, eps, f);
            if ~isempty(obj.X_last)
                del_x = X - obj.X_last;
                del_g = g - obj.g_last;
                
                obj.H = obj.H + (1 + (del_g'*obj.H*del_g)/(del_g'*del_x))*(del_x*del_x')/(del_x'*del_g)...
                              - (obj.H*del_g*del_x' + (obj.H*del_g*del_x')')/(del_g'*del_x);
            end
            
            p = -obj.H*g;
            obj.X_last = X;
            obj.g_last = g;
        end
        
        function Xk = test(obj)
            Xk = isempty(obj.X_last);
        end
    end
end

