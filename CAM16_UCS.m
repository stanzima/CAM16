% Author: Tanzima Habib, NTNU, 2025
% Github : https://github.com/stanzima/CAM16

function [J_ucs, a_ucs, b_ucs, h_ucs] = CAM16_UCS(J, M, h)
    % CAM16-UCS Uniform Color Space based on appendix B
    % in paper https://www.researchgate.net/publication/318152296_Comprehensive_color_solutions_CAM16_CAT16_and_CAM16-UCS cite:`Li2017`
    % Inputs:
    %   J: Lightness (from CAM16)
    %   M: Colorfulness (from CAM16)
    %   h: Hue angle (in degrees from CAM16)
    % Outputs:
    %   J_ucs: Lightness in CAM16-UCS
    %   a_ucs: a* coordinate in CAM16-UCS
    %   b_ucs: b* coordinate in CAM16-UCS
    %   h_ucs: Hue angle in CAM16-UCS (same as h)

    % Step 1: Transform lightness J to J_ucs
    J_ucs = (1.7 * J) / (1 + 0.007 * J);

    % Step 2: Transform colorfulness M to M_ucs
    M_ucs = log(1 + 0.0228 * M) / 0.0228;

    % Step 3: Compute a_ucs and b_ucs coordinates
    a_ucs = M_ucs * cosd(h); % a* coordinate
    b_ucs = M_ucs * sind(h); % b* coordinate

    % Step 4: Hue angle h_ucs remains the same
    h_ucs = h;
end
