
function out = SENSE_3D(input,sens)
    R=1;
    [Nx,Ny,Nz,~] = size(input);
    out = zeros(Nx,Ny,Nz);
    % loop over the top-1/R of the image
    for x = 1:Nx/R
        x_idx = x:Nx/R:Nx;
        % loop over the entire left-right extent
        for y = 1:Ny
            for z = 1:Nz
                % pick out the sub-problem sensitivities
                S = transpose(reshape(sens(x_idx,y,z,:),1,[]));
                % solve the sub-problem in the least-squares sense
                out(x_idx,y,z) = pinv(S)*reshape(input(x,y,z,:),[],1);
            end
        end
    end
end