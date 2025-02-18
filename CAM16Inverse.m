% Author: Tanzima Habib, NTNU, 2025
% Github : https://github.com/stanzima/CAM16


function XYZ = CAM16Inverse(J, C, h, XYZ_w, L_A, Y_b, surround)
    % CAM16 Color Appearance Model (Inverse Model) based on appendix A
    % in paper https://www.researchgate.net/publication/318152296_Comprehensive_color_solutions_CAM16_CAT16_and_CAM16-UCS cite:`Li2017`
    % Inputs:
    %   J: Lightness
    %   C: Chroma
    %   h: Hue angle (degrees)
    %   XYZ_w: 3x1 vector of the reference white in XYZ
    %   L_A: Adapting luminance in cd/m²
    %   Y_b: Background luminance factor (e.g., 20 for 20%)
    %   surround: Surround condition ('average', 'dim', 'dark')
    % Output:
    %   XYZ: 3x1 vector of the test color in XYZ 

    % Variables
    % RGB are different cone responses

    % Surround parameters
    switch surround
        case 'average'
            F = 1.0;
            c = 0.69;
            N_c = 1.0;
        case 'dim'
            F = 0.9;
            c = 0.59;
            N_c = 0.9;
        case 'dark'
            F = 0.8;
            c = 0.525;
            N_c = 0.8;
        otherwise
            error('Invalid surround. Use average, dim, or dark.');
    end

    % CAT16 transformation matrix
    M_CAT16 = [0.401288, 0.650173, -0.051461;
              -0.250268, 1.204414, 0.045854;
              -0.002079, 0.048952, 0.953127];

    % Inverse of CAT16 matrix
    M_CAT16_inv = [1.86206786, -1.01125463, 0.14918677;
                   0.38752654, 0.62144744, -0.00897398;
                  -0.01584150, -0.03412294, 1.04996444];

    % Step 0: Calculate viewing parameters
    RGB_w = M_CAT16 * XYZ_w;
    D = F * (1 - (1/3.6) * exp((-L_A - 42)/92));
    D = max(min(D, 1.0), 0.0); % If D is greater than one or less than zero, set it to one or zero, respectively

    D_R = D * (XYZ_w(2) / RGB_w(1)) + 1 - D;
    D_G = D * (XYZ_w(2) / RGB_w(2)) + 1 - D;
    D_B = D * (XYZ_w(2) / RGB_w(3)) + 1 - D;

    k = 1 / (5 * L_A + 1);
    F_L = 0.2 * k^4 * (5 * L_A) + 0.1 * (1 - k^4)^2 * (5 * L_A)^(1/3);

    n = Y_b / XYZ_w(2);
    z = 1.48 + sqrt(n);
    N_bb = 0.725 * (1 / n)^0.2;
    N_cb = N_bb;

    RGB_wc = [D_R * RGB_w(1); D_G * RGB_w(2); D_B * RGB_w(3)];
    RGB_aw = (400 * (F_L * RGB_wc / 100).^0.42) ./ ((F_L * RGB_wc / 100).^0.42 + 27.13) + 0.1;

    A_w = (2 * RGB_aw(1) + RGB_aw(2) + 0.05 * RGB_aw(3) - 0.305) * N_bb;

    % Step 1: Calculate t, e_t, A, p1, p2, p3
    t = (C / ((J / 100)^0.5 * (1.64 - 0.29^n)^0.73))^(1 / 0.9);
    e_t = (cos(h*pi/180 + 2) + 3.8) / 4;
    A = A_w * (J / 100)^(1 / (c * z));

    p1 = (50000 / 13) * N_c * N_cb * e_t / t;
    p2 = A / N_bb + 0.305;
    p3 = 21 / 20;

    % Step 2: Calculate a and b
    if t == 0
        a = 0;
        b = 0;
    else
        if abs(sind(h)) >= abs(cosd(h))
            p4 = p1 / sind(h);
            b = p2 * (2 + p3) * (460 / 1403) / (p4 + (2 + p3) * (220 / 1403) * (cosd(h) / sind(h)) - (27 / 1403) + p3 * (6300 / 1403));
            a = b * (cosd(h) / sind(h));
        else
            p5 = p1 / cosd(h);
            a = p2 * (2 + p3) * (460 / 1403) / (p5 + (2 + p3) * (220 / 1403) - (27 / 1403) * (sind(h) / cosd(h)) + p3 * (6300 / 1403));
            b = a * (sind(h) / cosd(h));
        end
    end

    % Step 3: Calculate R_a, G_a, B_a the post-adaptation cone responses
    R_a = (460 / 1403) * p2 + (451 / 1403) * a + (288 / 1403) * b;
    G_a = (460 / 1403) * p2 - (891 / 1403) * a - (261 / 1403) * b;
    B_a = (460 / 1403) * p2 - (220 / 1403) * a - (6300 / 1403) * b;

    % Step 4: Calculate R_c, G_c, B_c (cone response in the corresponding cone space due to the color adaptation of the illuminant) from post adaptation cone responses by undoing the dynamic range compression
    R_c = sign(R_a - 0.1) * (100 / F_L) * ((27.13 * abs(R_a - 0.1)) / (400 - abs(R_a - 0.1)))^(1 / 0.42);
    G_c = sign(G_a - 0.1) * (100 / F_L) * ((27.13 * abs(G_a - 0.1)) / (400 - abs(G_a - 0.1)))^(1 / 0.42);
    B_c = sign(B_a - 0.1) * (100 / F_L) * ((27.13 * abs(B_a - 0.1)) / (400 - abs(B_a - 0.1)))^(1 / 0.42);

    % Step 5: Calculate R, G, B cone responses from the corresponding cone responses by discounting the various luminance levels and surround conditions included in D
    RGB = [R_c / D_R; G_c / D_G; B_c / D_B];

    % Step 6: Calculate X, Y, Z
    XYZ = M_CAT16_inv * RGB;
end
