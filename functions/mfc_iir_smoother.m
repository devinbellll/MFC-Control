function x_filt = mfc_iir_smoother(x_raw, x_filt_km1, x_filt_km2, window)
%MFC_IIR_SMOOTHER Second-order recursive smoother shared by the MFC blocks.
%
%   x_filt = MFC_IIR_SMOOTHER(x_raw, x_filt_km1, x_filt_km2, window)
%
%   Unity-DC-gain, critically damped, second-order IIR low-pass:
%
%       x_filt[k] = ( x_raw[k] + (2W^2 + 2W) x_filt[k-1] - W^2 x_filt[k-2] )
%                   / (W^2 + 2W + 1)
%
%   where W = window is dimensionless (memory length in samples). W = 0 is
%   an exact pass-through. The filter has a repeated real pole at
%   z = W/(W+1), i.e. a time constant of roughly W samples.
%
%   Used for:
%     * the reference trajectory filter of the MFC controller, and
%     * the numerator/denominator smoothing of the algebraic F estimators.
%
%   Inputs
%     x_raw       : current raw sample
%     x_filt_km1  : previous filtered output   x_filt[k-1]
%     x_filt_km2  : filtered output before it  x_filt[k-2]
%     window      : W >= 0, memory in samples
%
%   See also MFC_FHAT_ALGEBRAIC_FIRST_ORDER, MFC_FHAT_ALGEBRAIC_SECOND_ORDER,
%   MFC_SISO_STEP.

x_filt = (x_raw + (2*window^2 + 2*window)*x_filt_km1 - window^2*x_filt_km2) ...
         / (window^2 + 2*window + 1);
end
