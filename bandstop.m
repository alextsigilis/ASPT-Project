% Function Signature:
% ----------------------------------------------
% y = bandstop(x, fs, m, s)
% ==============================================
%
% Function Description:
% ----------------------------------------------
% This function uses a Gaussian Filter to
% suppress an unwanted frequency band in
% some input signal x(t).
%
% In theory, the suppression of the unwanted
% frequencies could be achieved using an
% ideal band-cut filter with a rectangular
% transfer function in the frequency domain.
% However, such techniques are rarely used
% in practice, because they introduce ringing
% artifacts in the time-domain related to
% the Gibbs phenomenon.
%
% Therefore, this suppression is achieved by
% multiplying the spectrum of the input signal
% with an inverted Gaussian curve centered
% around the unwanted frequency m with a
% standard deviation of s.
%
% The transfer function of the filter has the
% following form:
% H(f) = 1 - exp(-0.5*((f-m)/s)^2)
%
% The spectrum of the input signal is obtained
% by applying the FFT algorithm on x(t):
% X(f) = FFT{x(t)}
%
% Filtering is carried out in the frequency
% domain through point-wise multiplication of
% X(f) and H(f):
% Y(f) = H(f) .* X(f)
%
% Finally, we return in the time domain by
% applying the IFFT algorithm on Y(f):
% y(t) = IFFT{Y(f)}
%
% For the purposes of our project, this function
% will be used to remove powerline interference
% from EEG, EOG and ECG signals. This type of
% interference manifests itself in the frequency
% domain as a spike in the 50Hz or 60Hz range
% depending on the frequency of the powerlines.
% ==============================================
%
% Arguments List:
% ----------------------------------------------
% x: (vector/1D array) the input signal. In our
% case an EEG/EOG/ECG signal contaminated with
% powerline interference.
%
% fs: (float/double) the sampling frequency
% of x expressed in Hertz (samples per second).
% In our case this is 256Hz.
%
% m: (float/double) the frequency component
% which is to be removed from x expressed in
% Hertz. In our case this is 50Hz.
%
% s: (float/double) the standard deviation of
% the gaussian bell expressed in Hertz
% ==============================================
%
% Return List:
% ----------------------------------------------
%
% y: (vector/1D-array)the result of applying
% the specified gaussian filter in the input
% sequence x.
% ==============================================
%
% Warnings and Limitations:
% 1) If x is a row vector, then y is
% a row vector as well.
%
% 2) If x is a column vector, then y is
% a column vector as well.
%
% 3) If x is a 2D array, then every column
% is treated as a separate signal. The same
% filter is applied on every column of the
% array. The return variable y will be a
% 2D array of the same dimensions as x.
%
% 4) The length the input signal(s)
% must be an even number.
%
% 5) The value of the cut-off frequency must
% be less than half of the sampling frequency.

function y = bandstop (x, fs, m, s)
  % 1) Transition to
  % the frequency domain
  x = x(:);
  X = fft(x);                           % apply the FFT algorithm
  X = fftshift(X);                      % shift zero frequency in the middle
  N = length(X);                        % length of FFT output

  % 2) Construct the transfer
  % function of the gaussian filter
  % and apply frequency domain
  % filtering
  k = N/2+1;                            % index of zero frequency element

  f = (fs/N)*[0:1:N/2-1]; f = f(:);     % positive frequencies
  H = 1 - exp(-0.5*((f-m)/s).^2);       % transfer function for positive axis
  X(k:N) = H .* X(k:N);                 % apply filter on positive axis

  f = (fs/N)*[-N/2:1:-1]; f = f(:);     % negative frequencies
  H = 1 - exp(-0.5*((f+m)/s).^2);       % transfer function for negative axis
  X(1:k-1) = H .* X(1:k-1);             % apply filter on negative axis

  % 4) Transition back
  % to the time domain
  X = ifftshift(X);                     % reverse the effect of fftshift()
  y = ifft(X);                          % apply the IFFT algorithm
  y = real(y);                          % discard rounding errors
end
