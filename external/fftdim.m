% https://htmlpreview.github.io/?https://github.com/mchiew/partial-fourier-tutorial/blob/main/PartialFourier.html
function x = fftdim(x,dims)
    for i = dims
        x = fftshift(fft(ifftshift(x,i),[],i),i);
    end
end
