function [Y] = unstagger(X,axis)
%UNSTAGGER Summary of this function goes here
%   Detailed explanation goes here
dim = sum(size(X)~=1);
if dim == 1
    Y = (X(1:end-1) + X(2:end))/2;
elseif dim == 2
    if axis == 0
        Y = (X(1:end-1,:) + X(2:end,:))/2;
    elseif axis == 1
        Y = (X(:,1:end-1) + X(:,2:end))/2;
    end
elseif dim == 3
    if axis == 0
        Y = (X(1:end-1,:,:) + X(2:end,:,:))/2;
    elseif axis == 1
        Y = (X(:,1:end-1,:) + X(:,2:end,:))/2;
    elseif axis == 2
        Y = (X(:,:,1:end-1) + X(:,:,2:end))/2;
    end
elseif dim == 4
    if axis == 0
        Y = (X(1:end-1,:,:,:) + X(2:end,:,:,:))/2;
    elseif axis == 1
        Y = (X(:,1:end-1,:,:) + X(:,2:end,:,:))/2;
    elseif axis == 2
        Y = (X(:,:,1:end-1,:) + X(:,:,2:end,:))/2;
    elseif axis == 3
        Y = (X(:,:,:,1:end-1) + X(:,:,:,2:end))/2;
    end
end
%}
end

