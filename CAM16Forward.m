% Author: Tanzima Habib, NTNU, 2025
% Github : https://github.com/stanzima/CAM16

function [J, Q, M, s, h, C, H, Hc] = CAM16Forward(XYZ, XYZ_w, L_A, Y_b, surround)
    % CAM16 Colour Appearance Model (Forward Model) based on appendix A
    % in paper https://www.researchgate.net/publication/318152296_Comprehensive_color_solutions_CAM16_CAT16_and_CAM16-UCS cite:`Li2017`
    % Inputs:
    %   XYZ: 3x1 vector of the test color in XYZ (scale [0, 100])
    %   XYZ_w: 3x1 vector of the reference white in XYZ
    %   L_A: Adapting luminance in cd/m²
    %   Y_b: Background luminance factor (e.g., 20 for 20%)
    %   surround: Surround condition ('average', 'dim', 'dark')
    % Outputs:
    %   J: Lightness
    %   Q: Brightness
    %   M: Colorfulness
    %   s: Saturation
    %   h: Hue angle (degrees)
    %   C: Chroma
    %   H: Hue composition
    %   Hc:Hue composition report e.g.: Hc ='G76B24' means 76% Green and  24% Blue

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

    % Step 0: Calculate all values/parameters independent of the input sample
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

    RGB_wc = zeros(1,3);
    RGB_wc(1) = D_R * RGB_w(1);
    RGB_wc(2) = D_G * RGB_w(2); 
    RGB_wc(3) = D_B * RGB_w(3);

    RGB_aw = (400 * ((F_L * RGB_wc / 100).^0.42 ./ ((F_L * RGB_wc / 100).^0.42 + 27.13))) + 0.1;

    A_w = (2 * RGB_aw(1) + RGB_aw(2) + (RGB_aw(3)/20) - 0.305) * N_bb;

    % Step 1: Calculate cone responses
    RGB = M_CAT16 * XYZ;

    % Step 2: Complete the color adaptation
    RGB_c = [D_R * RGB(1); D_G * RGB(2); D_B * RGB(3)];

    % Step 3: Calculate post-adaptation cone responses
    RGB_a = (400 * ((F_L * RGB_c / 100).^0.42 ./ ((F_L * RGB_c / 100).^0.42 + 27.13))) + 0.1;
    RGB_a(RGB_c<0) = -400 * ((-1*F_L * RGB_c(RGB_c<0) / 100).^0.42 ./ ((-1*F_L * RGB_c(RGB_c<0) / 100).^0.42 + 27.13)) + 0.1;

    % Step 4: Calculate a, b, and hue angle h
    a = RGB_a(1) - 12/11 * RGB_a(2) + 1/11 * RGB_a(3);
    b = (1/9) * (RGB_a(1) + RGB_a(2) - 2 * RGB_a(3));
    h = atan2d(b, a);
    if h < 0
        h = h + 360;
    end
    tableA_h = [20.14; 90; 164.25; 237.53; 380.14]; 
    tableA_e = [0.8; 0.7; 1; 1.2; 0.8];
    tableA_H = [0; 100; 200; 300; 400];
    hue_text = ['R', 'Y', 'G', 'B', 'R'];
    

    if h < tableA_h(1)  %Set unique hue data from Table A2
        i = 3;
        h = h + 360;
        
    elseif h >= tableA_h(1) && h<tableA_h(2)
        i = 1;
    elseif h >= tableA_h(2) && h<tableA_h(3)
        i = 2; 
    elseif h >= tableA_h(3) && h<tableA_h(4)
        i = 3;
    elseif h >= tableA_h(4) && h<tableA_h(5)
        i = 4;   
    end


    % Step 5: Calculate eccentricity and hue composition
    e_t = (cos(h*pi/180 + 2) + 3.8) / 4;
    H = tableA_H(i) + (100*(h-tableA_h(i))/tableA_e(i))/((h-tableA_h(i))/tableA_e(i) + (tableA_h(i+1)-h)/tableA_e(i+1)); %Hue Quadrature H
    
    
    Hc = strcat(hue_text(i),num2str(round(H-tableA_H(i))),hue_text(i+1),num2str(round(tableA_H(i+1)-H))); %Hue composition report e.g.: Hc ='G76B24' means 76% Green and  24% Blue 

    % Step 6: Calculate achromatic response A
    A = (2 * RGB_a(1) + RGB_a(2) + (RGB_a(3)/20) - 0.305) * N_bb;

    % Step 7: Calculate correlate of lightness J
    J = 100 * (A / A_w)^(c * z);

    % Step 8: Calculate correlate of brightness Q
    Q = (4 / c) * sqrt(J / 100) * (A_w + 4) * F_L^0.25;

    % Step 9: Calculate correlates of chroma C, colorfulness M, and saturation s
    t = ((50000 / 13 * N_c * N_cb * e_t) * sqrt(a^2 + b^2)) / (RGB_a(1) + RGB_a(2) + 21/20 * RGB_a(3));
    C = t^0.9 * sqrt(J / 100) * (1.64 - 0.29^n)^0.73;
    M = C *(F_L)^0.25;
    s = 100 * sqrt(M / Q);
end
