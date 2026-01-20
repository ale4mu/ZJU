% Example 4.1
function [U_L, U_R] = WENO_gPC(U_hat, Phi_val, gamma)
    % U_hat：单元平均值
    % 用单元内的gauss-labtuo节点重构并做拉格朗日插值
    %然后计算左边界处的右极限和左极限
    [nVar, nModes, Nx] = size(U_hat);
    
    % 周期性边界
    idx_ext = [Nx-2, Nx-1, Nx, 1:Nx, 1, 2, 3];
    U_ext = U_hat(:, :, idx_ext);
    
    U_L = zeros(nVar, nModes, Nx+1);
    U_R = zeros(nVar, nModes, Nx+1);
    
    for i = 1:Nx
        idx = i + 3; 
        for k = 1:nModes
            for v = 1:nVar
                % [i-2 ... i+2]
                stencil = squeeze(U_ext(v, k, idx-2 : idx+2));
                % 正通量：单元i在右边界的值
                U_L(v,k, i+1) = WENO_pos(stencil); 
                % 负通量：单元i+1在左边界的值
                stencil_neg = squeeze(U_ext(v, k, idx-1 : idx+3));
                U_R(v, k, i+1) = WENO_neg(stencil_neg);
            end
        end
    end

    % 周期性
    U_L(:,:,1) = U_L(:,:,Nx+1);
    U_R(:,:,1) = U_R(:,:,Nx+1);
end