function [Normalization] = Norm(Nr)
        

        Str = load('Xmn.mat');
        Xmn = Str.Xmn;

        for i = 1:length(Nr)
            disp(Nr(i));
            mode = Xmn(Nr(i)).mode;
            m = Xmn(Nr(i)).m;
            xmn = Xmn(Nr(i)).xmn;
            
            if m == 0
                deltan = 1;
            else
                deltan = 0;
            end
            
            if mode == "TE"
                
                Normalization(i, :) = abs(sqrt(2/pi) ./ (besselj(m, xmn) .* sqrt(xmn.^2 - m.^2) .* sqrt(1 + deltan)));
                
            elseif mode == "TM"
                
                Normalization(i, :) = abs(sqrt(2/pi) ./ (besselj_der(m, xmn) .* xmn .* sqrt(1 + deltan)));
                
            end
            
        end

end

