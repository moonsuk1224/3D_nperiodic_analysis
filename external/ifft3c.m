function res = ifft3c(x)
    S = size(x);
    if length(S) < 3
        S(3) = 1;  % Handle 2D case
    end
    fctr = sqrt(S(1)*S(2)*S(3));
    
    % Handle different input dimensions
    if ndims(x) == 2
        res = ifft2c(x);  % Use your existing 2D function
    elseif ndims(x) == 3
        res = fctr * fftshift(ifftn(ifftshift(x)));
    else
        % For 4D+ arrays, apply 3D IFFT to first 3 dimensions
        orig_size = size(x);
        x_reshaped = reshape(x, S(1), S(2), S(3), []);
        res = zeros(size(x_reshaped), 'like', x);
        for i = 1:size(x_reshaped, 4)
            res(:,:,:,i) = fctr * fftshift(ifftn(ifftshift(x_reshaped(:,:,:,i))));
        end
        res = reshape(res, orig_size);
    end
end