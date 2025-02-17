% Author: Tanzima Habib, NTNU, 2025
% Github : https://github.com/stanzima/CAM16

function [DECAM16UCS, CIELAB] = CAM16_UCS_Colour_Difference(J_ucs1, a_ucs1, b_ucs1, J_ucs2, a_ucs2, b_ucs2)
    % Compute Euclidean distance (colour difference) in CAM16-UCS based on appendix B
    % in paper https://www.researchgate.net/publication/318152296_Comprehensive_color_solutions_CAM16_CAT16_and_CAM16-UCS cite:`Li2017`
    % Inputs:
    %   J_ucs1, a_ucs1, b_ucs1: CAM16-UCS coordinates of colour sample 1
    %   J_ucs2, a_ucs2, b_ucs2: CAM16-UCS coordinates of colour sample2
    % Output:
    %   DECAM16UCS(Delta E CAM16 UCS): CAM16-UCS Colour difference i.e. euclidean distance in CAM16-UCS

    % Step 1: Compute differences in J_ucs, a_ucs, and b_ucs
    delta_J_ucs = J_ucs2 - J_ucs1;
    delta_a_ucs = a_ucs2 - a_ucs1;
    delta_b_ucs = b_ucs2 - b_ucs1;

    % Step 2: Compute Euclidean distance
    DECAM16UCS = sqrt(delta_J_ucs^2 + delta_a_ucs^2 + delta_b_ucs^2);

    % Step 2: Compute CIELAB
    CIELAB = 1.41*(DECAM16UCS)^0.63;

end