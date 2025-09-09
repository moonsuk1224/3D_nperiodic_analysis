% https://htmlpreview.github.io/?https://github.com/mchiew/partial-fourier-tutorial/blob/main/PartialFourier.html
function x = ifftdim(x,dims)
    for i = dims
        x = fftshift(ifft(ifftshift(x,i),[],i),i);
    end
end