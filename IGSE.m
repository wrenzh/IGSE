% IGSE core loss calculation
%
% Author: Wen Zhang (wren.zh@gmail.com) 
% Version: 0.2 
% License: MIT 
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
% WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
% COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% This script implements the improved general Steinmetz equation core loss calculation method. For more information on
% the techincal details, please see DOI: 10.1109/CIPE.2002.1196712.
%
% Inputs:
%     time : evenly sampled time vector in seconds
%     flux : flux density vector in seconds
%     a    : Steinmetz coefficient alpha
%     b    : Steinmetz coefficient beta
%     k    : Steinmetz coefficient k
%     tol  : flux periocity check tolerance (optional, 0.1)
%
% Output:
%     loss : core loss in W/m^3
%
% The optional argument "tol" is the tolerance used when checking flux periodicity. Ideally, the input should be a
% single cycle of period signal and flux(1) == flux(end). However, if the input is from simulation, it is more often
% than not that this is not true. This argument suggest how different is tolerable, with 0 being the most strict and 1
% being the most loose.
%
% The best execution speed is about 60x-120x faster than previous implementation. For best execution speed, provide a
% over-sampled single cycle of time and flux signal with even time spacing. This is usually true for simulation results
% coming from MATLAB/Simulink. For other simulators, consider resample the signal in MATLAB first.
%
% Caveat: the implementation here does NOT perform interpolation at minor loop boundaries for the sake of speed. For
% most cases where the signal is sufficiently sampled like in simulation, there isn't any issue. However, if the input
% is not over-sampled, the error may be significant.

function loss = IGSE(time, flux, a, b, k, tol)
    arguments
        time (1, :) double
        flux (1, :) double
        a    (1, 1) double {mustBePositive}
        b    (1, 1) double {mustBePositive}
        k    (1, 1) double {mustBePositive}
        tol  (1, 1) double {mustBePositive} = 0.1
    end
    
    % time vector should be in increasing order
    if any(diff(time) < 0)
        error('Invalid input argument. Time vector must be successive.');
    end
    
    % time vector should be evenly sampled
    if any(diff(diff(time)) > 0.01 * range(time))
        error('Invalid input argument. Time vector should be evenly spaced. Consider resample inputs.');
    end
    
    % time and flux vector must have the same length
    if numel(time) ~= numel(flux)
        error('Invalid input argument. Time and flux vectors must have same number of elements.');
    end
    
    % flux vector should be periodic so flux(1) == flux(end)
    if abs(flux(1) - flux(end)) > tol * range(flux)
        error('Invalid input argument. Flux vector must be periodic.');
    end
    
    % Reshaping input vectors to column vectors
    time = time(:);
    flux = flux(:);
    
    % Reorganize so the flux vector starts from minimum
    [~, idx] = min(flux);
    flux = [flux(idx : end); flux(2 : idx)];
    time = [time(idx : end); time(2 : idx) + time(end) - time(1)];
    
    % Average the energy loss over the entire time period
    loss = 1000 * sumMinorLoops(time, flux, a, b, k) / range(time);
end


% Sum all minor loop energy loss together
%
% Inputs:
%    t : time vector
%    B : flux vector
%    a : Steinmetz coefficient a
%    b : Steinmetz coefficient b
%    k : Steinmetz coefficient k
%
% Output:
%    e : core energy loss in mW/cm^3
function e = sumMinorLoops(t, B, a, b, k)
    ki = modifiedCoefficient(a, b, k);
    marks = markMinorLoops(B);
    dt = t(2) - t(1);
    
    e = 0;
    for m = 0 : max(marks)
        Bm = B(marks == m);
        e = e + ki * range(Bm) ^ (b - a) * sum(abs(diff(Bm) / dt).^a * dt);
    end
end


% Identify the start and end indices of minor loop by marking a linked array
%
% Input:
%    B is the orignal flux vector
%
% Output:
%    marks is an array with the same number of elements of B
% 
% The array is initialized to be all zeros. At the first outmost minor loop start and end indices, the value of the
% array is marked as one. Any immediate nested minor loop or subsequent minor loop will be marked as two, and so on.
%
% Note the marks are stored in uint32 so there is a theoretical depth limit to minor loop nesting.
% Also note the flux vector should be of sufficient detail. Choppy waveforms with nothing but ups and downs
% will result in meaningless calculation.
function marks = markMinorLoops(B)
    marks = zeros(1, numel(B), 'uint32');
    d = diff(B);
    
    loopNum = 1;
    region = [1, numel(B)];
    
    while true
        [sp, ep] = outerMostMinorLoop(B, d, region);
        
        % No more minor loops in this region, move on
        if sp == 0 && ep == 0
            % Search ends at the end of vector
            if region(2) == numel(B)
                break;
            % Continue searching in the immediate parent region after the current region
            else
                region(1) = region(2) + 1;
                region(2) = find(marks == marks(region(2)+1), 1, 'last');
            end
        % Minor loops found in this region, check for nested minor loops
        else
            % Mark and increment 
            marks(sp : ep) = loopNum;
            
            loopNum = loopNum + 1;
            
            % Check for integer overflow
            if loopNum == intmax('uint32')
                error('Runtime error. Too many minor loops.');
            end
            
            region = [sp+1, ep];
        end
    end
end

% Find the start and end indices of the first outermost minor loop.
%
% Inputs:
%    B is the original flux vector
%    d is the diff(B) and is passed instead of calculated in-place to save calculation
%    region is the region to inspect and is a two-element vector [startIndex, endIndex]
%
% Outputs:
%    [sp, ep] are the start and end points of the minor loop. If there is no minor loop, sp == ep == 0.
%
% Note when B(sp) is local maximum, B(ep) >= B(sp) and B(ep-1) < B(sp).
% And when B(sp) is local minimum, B(ep) <= B(sp) and B(ep-1) > B(sp).
function [sp, ep] = outerMostMinorLoop(B, d, region)
    sp = 0;
    ep = 0;
    
    for i = region(1)+1 : region(2)-1
        % Local maxima
        if d(i-1) > 0 && d(i) < 0
            len = find(B(i+1 : end) >= B(i), 1);
            if ~isempty(len) && i + len <= region(2)
                sp = i;
                ep = i + len;
                break;
            end
        end
        
        % Local minima
        if d(i-1) < 0 && d(i) > 0 && any(B(i+1 : end) <= B(i))
            len = find(B(i+1 : end) <= B(i), 1);
            if ~isempty(len) && i + len <= region(2)
                sp = i;
                ep = i + len;
                break;
            end
        end
    end
end


% Numerical integration for the modified Steinmetz coefficient ki
function ki = modifiedCoefficient(a, b, k)
    x = linspace(0, pi / 2);
    q = 4 * trapz(x, cos(x) .^ a);
    ki = k / (2^(b-1) * pi^(a-1) * q);
end