% UBIT Name   : ankitary
% Name        : Ankit Arya
% UB Person # : 5009-6802

function [] = Project()                                                    % MAIN FUNCTION %
    clc;
    I = imread('IMG_0162.JPG');
    J = rgb2gray(I);
    imwrite(J, 'Gray3.jpg');
    B = imresize(J, [1024 1024]);
    imwrite(B, 'Gray3Sub.jpg');
    x = double(B);
    y = uint8(zeros(size(B)));
    globe = 0;
    qjpg = [16 11 10 16  24  40  51  61
            12 12 14 19  26  58  60  55
            14 13 16 24  40  57  69  56
            14 17 22 29  51  87  80  62
            18 22 37 56  68 109 103  77
            24 35 55 64  81 104 113  92
            49 64 78 87 103 121 120 101
            72 92 95 98 112 100 103  99];
    final = zeros(1024, 1024);
    ct = 0;
    cg = 0;
    cs = 0;
    binary = [];
    sending = [];
    for r = 1:8:1024
        for c = 1:8:1024
            if ((c > 8) || (r > 8)) && numel(HH) ~= 0 
                cg = HH(1);
            else
                cg = 0; 
            end
            x = B(r:r+7, c:c+7)-128;
            V = dct2(x);
            V = round(V./qjpg).*qjpg;
            FF = zigizag(V);
            [HH temp] = run_length(FF); %#ok<*ASGLU,*NASGU>
            HH(1) = cg - temp;
            if HH(1) < 0
                HH(1) = -(HH(1));
            end
            bitstream = encoder(HH);
            for num = 1 : length(bitstream)
                tep = str2num(bitstream(num)); %#ok<ST2NM>
                binary = [binary, tep];
            end
            count = sizeof_bitstream(bitstream);
            cs = cs + count;   % For counting the number of bits in the bitstream for the whole image to caculate the compression ratio %
            if(r > 1016 && c > 1016)
                for ii = 1:3:size(binary,2)
                    fg = binary(ii:ii+2);
                    lbc_encoded = lbc_encoder(fg);
                    lbc_decoded = lbc_decoder(fg,lbc_encoded);
                    sending = [sending, lbc_encoded];
                end
            end
            FC = Inverse_runlength(HH);
            DE = Izigizag(FC);
            FD = (DE.*qjpg);
            final(r:r+7,c:c+7)=FD;
            ct = ct + sum(DE(:)~=0);
            y(r:r+7, c:c+7) = idct2(DE)+128;
       end
    end
    compression_ratio = cs/(1024*1024*8);
    fprintf('The Compression Ratio is : ');
    disp(compression_ratio);
    figure; imshow(B);
    figure; imshow(y);
    imwrite(y, 'Gray3Comp_Standard_JPEG.jpg');
    [row column] = size(B);
    mse = (double(B) - double(y)).^2;
    mm = sum(sum(mse))/(row * column);
    PSNR_Value = 10 * log10(255^2/mm);    
    fprintf('PSNR VALUE IS : ');
    disp(PSNR_Value);
end

function KK = zigizag(V)                                                   % FUNCTION TO COMPUTE ZIGZAG SCAN %

    KK =   [ V(1,1) V(1,2) V(2,1) V(3,1) V(2,2) V(1,3) V(1,4) V(2,3) V(3,2) V(4,1)...
             V(5,1) V(4,2) V(3,3) V(2,4) V(1,5) V(1,6) V(2,5) V(3,4) V(4,3) V(5,2)...
             V(6,1) V(7,1) V(6,2) V(5,3) V(4,4) V(3,5) V(2,6) V(1,7) V(1,8) V(2,7)...
             V(3,6) V(4,5) V(5,4) V(6,3) V(7,2) V(8,1) V(8,2) V(7,3) V(6,4) V(5,5)...
             V(4,6) V(3,7) V(2,8) V(3,8) V(4,7) V(5,6) V(6,5) V(7,4) V(8,3) V(8,4)...
             V(7,5) V(6,6) V(5,7) V(4,8) V(5,8) V(6,7) V(7,6) V(8,5) V(8,6) V(7,7)...
             V(6,8) V(7,8) V(8,7) V(8,8)];
end

function YY = Izigizag(FF)                                                 % FUNCTION TO COMPUTE INVERSE-ZIGZAG SCAN %

    YY =    [ FF(1) FF(2) FF(6) FF(7) FF(15) FF(16) FF(28) FF(29);
              FF(3) FF(5) FF(8) FF(14) FF(17) FF(27) FF(30) FF(43);
              FF(4) FF(9) FF(13) FF(18) FF(26) FF(31) FF(42) FF(44);
              FF(10) FF(12) FF(19) FF(25) FF(32) FF(41) FF(45) FF(54);
              FF(11) FF(20) FF(24) FF(33) FF(40) FF(46) FF(53) FF(55);
              FF(21) FF(23) FF(34) FF(39) FF(47) FF(52) FF(56) FF(61);
              FF(22) FF(35) FF(38) FF(48) FF(51) FF(57) FF(60) FF(62);
              FF(36) FF(37) FF(49) FF(50) FF(58) FF(59) FF(63) FF(64)];
end

function irle = Inverse_runlength(AA)                                      % FUNCTION TO REVERSE THE RUN LENGTH ENCODING %
    ms = size(AA);
    length = ms(2);
    if(numel(AA) == 0)
        irle = zeros(1,64);
    else
        irle = zeros(1,64);
        irle(1) = AA(1);
        j = 2;
        for i = 2:length
            if( rem(i,2) == 0 )
                while(j <= AA(i))
                    j = j + 1;
                end
            else
                irle(j) = AA(i);
            end
        end
    end
end

function [rle temp] = run_length(SD) %#ok<*INUSD>                          % FUNCTION TO DO RUN LENGTH ENCODING % 
    rle = [];
    if numel(SD) ~= 0
        temp = SD(1);
    else
        rle = 0;
        temp = 0;
    end
    ms = size(SD);
    length = ms(2);
    rl = 0;
    if SD(1) > 0 || SD(1) < 0
        vcc = SD(1);
        rle = [rle SD(1)];
    end
    for i=2:length
        if SD(i) == 0
            rl = rl + 1;
        else
            vc = SD(i);
            rle = [rle rl SD(i)]; %#ok<*AGROW>
            rl = 0;
        end
    end 
end

function bitstream = encoder(rle)                                          % FUNCTION TO PRODUCE THE BITSTREAM USING HUFFMAN CODING FUNCTIONS HARD-CODED BELOW %
    ms = size(rle);
    length = ms(2);
    bitstream = [];
    rl = 0;
    if(numel(rle) == 0)
        rle = zeros(1,64);
    end
    if rle(1) > 0 || rle(1) < 0
        vcc = rle(1);
        fgg = Baseline_DC(vcc);
        bvv = Default_DC(fgg);
        if vcc >= 0
            hp = dec2bin(vcc);
        else
            ab = abs(vcc);
            htp = dec2bin(ab);
            si = size(htp);
            sf = bitcmp(ab,si);
            hp = dec2bin(sf);
        end
        bvc = strcat(bvv,hp);
        bitstream = [bitstream bvc];
    end
    for i = 2 : length
        if rle(i) == 0
            rl = rl + 1;
        else
            vc = rle(i);
            fg = Baseline_AC(vc);
            bv = Huff_Code(rl,fg);
            if vc >= 0
                hpp = dec2bin(vc);
            else
                hpp = dec2bin(typecast(int8(vc),'uint8'));
            end
            bvb = strcat(bv,hpp);
            bitstream = [bitstream bvb];
            rl = 0;
        end
    end
    if rle == 0
        bitstream = '1010';
    else
        bitstream = [bitstream '1010'];
    end
end

function count = sizeof_bitstream(bitstream)                               % FUNCTION TO CALCULATE THE NUMBER OF BITS IN THE BITSTREAM %
    count = 0;
    len = length(bitstream);
    for i = 1:len
        %d = size(bitstream(i),2);
        count = count + 1;
    end
end

function bing = Huff_Code(Run,Run_Category) %#ok<*STOUT,*DEFNU>            % FUNCTION FOR GENERATING HUFFMAN BITSTREAM %

    if(Run == 0 && Run_Category == 1)
        bing = '00';
    elseif(Run == 0 && Run_Category == 2)
        bing = '01';
    elseif(Run == 0 && Run_Category == 3)
        bing = '100';
    elseif(Run == 0 && Run_Category == 4)
        bing = '1011';
    elseif(Run == 0 && Run_Category == 5)
        bing = '11010';
    elseif(Run == 0 && Run_Category == 6)
        bing = '111000';
    elseif(Run == 0 && Run_Category == 7)
        bing = '1111000';
    elseif(Run == 0 && Run_Category == 8)
        bing = '1111110110';
    elseif(Run == 0 && Run_Category == 9)
        bing = '1111111110000010';
    elseif(Run == 0 && Run_Category == 10)
        bing = '1111111110000011';
    elseif(Run == 1 && Run_Category == 1)
        bing = '1100';
    elseif(Run == 1 && Run_Category == 2)
        bing = '11011';
    elseif(Run == 1 && Run_Category == 3)
        bing = '1111001';
    elseif(Run == 1 && Run_Category == 4)
        bing = '111110110';
    elseif(Run == 1 && Run_Category == 5)
        bing = '11111110110';
    elseif(Run == 1 && Run_Category == 6)
        bing = '1111111110000100';
    elseif(Run == 1 && Run_Category == 7)
        bing = '1111111110000101';
    elseif(Run == 1 && Run_Category == 8)
        bing = '1111111110000110';
    elseif(Run == 1 && Run_Category == 9)
        bing = '1111111110000111';
    elseif(Run == 1 && Run_Category == 10)
        bing = '1111111110001000';
    elseif(Run == 2 && Run_Category == 1)
        bing = '11100';
    elseif(Run == 2 && Run_Category == 2)
        bing = '11111001';
    elseif(Run == 2 && Run_Category == 3)
        bing = '1111110111';
    elseif(Run == 2 && Run_Category == 4)
        bing = '111111110100';
    elseif(Run == 2 && Run_Category == 5)
        bing = '1111111110001001';
    elseif(Run == 2 && Run_Category == 6)
        bing = '1111111110001010';
    elseif(Run == 2 && Run_Category == 7)
        bing = '1111111110001011';
    elseif(Run == 2 && Run_Category == 8)
        bing = '1111111110001100';
    elseif(Run == 2 && Run_Category == 9)
        bing = '1111111110001101';
    elseif(Run == 2 && Run_Category == 10)
        bing = '1111111110001110';
    elseif(Run == 3 && Run_Category == 1)
        bing = '111010';
    elseif(Run == 3 && Run_Category == 2)
        bing = '111110111';
    elseif(Run == 3 && Run_Category == 3)
        bing = '111111110101';
    elseif(Run == 3 && Run_Category == 4)
        bing = '1111111110001111';
    elseif(Run == 3 && Run_Category == 5)
        bing = '1111111110010000';
    elseif(Run == 3 && Run_Category == 6)
        bing = '1111111110010001';
    elseif(Run == 3 && Run_Category == 7)
        bing = '1111111110010010';
    elseif(Run == 3 && Run_Category == 8)
        bing = '1111111110010011';
    elseif(Run == 3 && Run_Category == 9)
        bing = '1111111110010100';
    elseif(Run == 3 && Run_Category == 10)
        bing = '1111111110010101';
    elseif(Run == 4 && Run_Category == 1)
        bing = '111011';
    elseif(Run == 4 && Run_Category == 2)
        bing = '1111111000';
    elseif(Run == 4 && Run_Category == 3)
        bing = '1111111110010110';
    elseif(Run == 4 && Run_Category == 4)
        bing = '1111111110010111';
    elseif(Run == 4 && Run_Category == 5)
        bing = '1111111110011000';
    elseif(Run == 4 && Run_Category == 6)
        bing = '1111111110011001';
    elseif(Run == 4 && Run_Category == 7)
        bing = '1111111110011010';
    elseif(Run == 4 && Run_Category == 8)
        bing = '1111111110011011';
    elseif(Run == 4 && Run_Category == 9)
        bing = '1111111110011100';
    elseif(Run == 4 && Run_Category == 10)
        bing = '1111111110011101';
    elseif(Run == 5 && Run_Category == 1)
        bing = '1111010';
    elseif(Run == 5 && Run_Category == 2)
        bing = '11111110111';
    elseif(Run == 5 && Run_Category == 3)
        bing = '1111111110011110';
    elseif(Run == 5 && Run_Category == 4)
        bing = '1111111110011111';
    elseif(Run == 5 && Run_Category == 5)
        bing = '1111111110100000';
    elseif(Run == 5 && Run_Category == 6)
        bing = '1111111110100001';
    elseif(Run == 5 && Run_Category == 7)
        bing = '1111111110100010';
    elseif(Run == 5 && Run_Category == 8)
        bing = '1111111110100011';
    elseif(Run == 5 && Run_Category == 9)
        bing = '1111111110100100';
    elseif(Run == 5 && Run_Category == 10)
        bing = '1111111110100101';
    elseif(Run == 6 && Run_Category == 1)
        bing = '1111011';
    elseif(Run == 6 && Run_Category == 2)
        bing = '111111110110';
    elseif(Run == 6 && Run_Category == 3)
        bing = '1111111110100110';
    elseif(Run == 6 && Run_Category == 4)
        bing = '1111111110100111';
    elseif(Run == 6 && Run_Category == 5)
        bing = '1111111110101000';
    elseif(Run == 6 && Run_Category == 6)
        bing = '1111111110101001';
    elseif(Run == 6 && Run_Category == 7)
        bing = '1111111110101010';
    elseif(Run == 6 && Run_Category == 8)
        bing = '1111111110101011';
    elseif(Run == 6 && Run_Category == 9)
        bing = '1111111110101100';
    elseif(Run == 6 && Run_Category == 10)
        bing = '1111111110101101';
    elseif(Run == 7 && Run_Category == 1)
        bing = '11111010';
    elseif(Run == 7 && Run_Category == 2)
        bing = '111111110111';
    elseif(Run == 7 && Run_Category == 3)
        bing = '1111111110101110';
    elseif(Run == 7 && Run_Category == 4)
        bing = '1111111110101111';
    elseif(Run == 7 && Run_Category == 5)
        bing = '1111111110110000';
    elseif(Run == 7 && Run_Category == 6)
        bing = '1111111110110001';
    elseif(Run == 7 && Run_Category == 7)
        bing = '1111111110110010';
    elseif(Run == 7 && Run_Category == 8)
        bing = '1111111110110011';
    elseif(Run == 7 && Run_Category == 9)
        bing = '1111111110110100';
    elseif(Run == 7 && Run_Category == 10)
        bing = '1111111110110101';
    elseif(Run == 8 && Run_Category == 1)
        bing = '111111000';
    elseif(Run == 8 && Run_Category == 2)
        bing = '111111111000000';
    elseif(Run == 8 && Run_Category == 3)
        bing = '1111111110110110';
    elseif(Run == 8 && Run_Category == 4)
        bing = '1111111110110111';
    elseif(Run == 8 && Run_Category == 5)
        bing = '1111111110111000';
    elseif(Run == 8 && Run_Category == 6)
        bing = '1111111110111001';
    elseif(Run == 8 && Run_Category == 7)
        bing = '1111111110111010';
    elseif(Run == 8 && Run_Category == 8)
        bing = '1111111110111011';
    elseif(Run == 8 && Run_Category == 9)
        bing = '1111111110111100';
    elseif(Run == 8 && Run_Category == 10)
        bing = '1111111110111101';
    elseif(Run == 9 && Run_Category == 1)
        bing = '111111001';
    elseif(Run == 9 && Run_Category == 2)
        bing = '1111111110111110';
    elseif(Run == 9 && Run_Category == 3)
        bing = '1111111110111111';
    elseif(Run == 9 && Run_Category == 4)
        bing = '1111111111000000';
    elseif(Run == 9 && Run_Category == 5)
        bing = '1111111111000001';
    elseif(Run == 9 && Run_Category == 6)
        bing = '1111111111000010';
    elseif(Run == 9 && Run_Category == 7)
        bing = '1111111111000011';
    elseif(Run == 9 && Run_Category == 8)
        bing = '1111111111000100';
    elseif(Run == 9 && Run_Category == 9)
        bing = '1111111111000101';
    elseif(Run == 9 && Run_Category == 10)
        bing = '1111111111000110';
    elseif(Run == 10 && Run_Category == 1)
        bing = '111111010';
    elseif(Run == 10 && Run_Category == 2)
        bing = '1111111111000111';
    elseif(Run == 10 && Run_Category == 3)
        bing = '1111111111001000';
    elseif(Run == 10 && Run_Category == 4)
        bing = '1111111111001001';
    elseif(Run == 10 && Run_Category == 5)
        bing = '1111111111001010';
    elseif(Run == 10 && Run_Category == 6)
        bing = '1111111111001011';
    elseif(Run == 10 && Run_Category == 7)
        bing = '1111111111001100';
    elseif(Run == 10 && Run_Category == 8)
        bing = '1111111111001101';
    elseif(Run == 10 && Run_Category == 9)
        bing = '1111111111001110';
    elseif(Run == 10 && Run_Category == 10)
        bing = '1111111111001111';
    else
        bing = '00';
    end
end

function [ac] = Baseline_AC(Coeff)                                         % FUNCTION TO COMPUTE AC COEFFICIENTS %

    if Coeff == -1 || Coeff == 1
        ac = 1;
    elseif Coeff == -3 || Coeff == -2 || Coeff == 2 || Coeff == 3
        ac = 2;
    elseif (Coeff >= -7 && Coeff <= -4) || (Coeff >= 4 && Coeff <= 7)
        ac = 3;
    elseif (Coeff >= -15 && Coeff <= -8) || (Coeff >= 8 && Coeff <= 15)
        ac = 4;
    elseif (Coeff >= -31 && Coeff <= -16) || (Coeff >= 16 && Coeff <= 31)
        ac = 5;
    elseif (Coeff >= -63 && Coeff <= -32) || (Coeff >= 32 && Coeff <= 63)
        ac = 6;
    elseif (Coeff >= -127 && Coeff <= -64) || (Coeff >= 64 && Coeff <= 127)
        ac = 7;
    elseif (Coeff >= -255 && Coeff <= -128) || (Coeff >= 128 && Coeff <= 255)
        ac = 8;
    elseif (Coeff >= -511 && Coeff <= -256) || (Coeff >= 256 && Coeff <= 511)
        ac = 9;
    elseif (Coeff >= -1023 && Coeff <= -512) || (Coeff >= 512 && Coeff <= 1023)
        ac = 10;
    else
        ac = 1;
    end
end

function [dc] = Baseline_DC(Coeff)                                         % FUNCTION TO CALCULATE DC COEFFICIENTS % 

    if Coeff == 0
        dc = 0;
    elseif Coeff == -1 || Coeff == 1
        dc = 1;
    elseif Coeff == -3 || Coeff == -2 || Coeff == 2 || Coeff == 3
        dc = 2;
    elseif (Coeff >= -7 && Coeff <= -4) || (Coeff >= 4 && Coeff <= 7)
        dc = 3;
    elseif (Coeff >= -15 && Coeff <= -8) || (Coeff >= 8 && Coeff <= 15)
        dc = 4;
    elseif (Coeff >= -31 && Coeff <= -16) || (Coeff >= 16 && Coeff <= 31)
        dc = 5;
    elseif (Coeff >= -63 && Coeff <= -32) || (Coeff >= 32 && Coeff <= 63)
        dc = 6;
    elseif (Coeff >= -127 && Coeff <= -64) || (Coeff >= 64 && Coeff <= 127)
        dc = 7;
    elseif (Coeff >= -255 && Coeff <= -128) || (Coeff >= 128 && Coeff <= 255)
        dc = 8;
    elseif (Coeff >= -511 && Coeff <= -256) || (Coeff >= 256 && Coeff <= 511)
        dc = 9;
    elseif (Coeff >= -1023 && Coeff <= -512) || (Coeff >= 512 && Coeff <= 1023)
        dc = 10;
    else
        dc = 0;
    end
end

function huffman_dc = Default_DC(dc_coeff)                                 % FUNCTION FOR DEFAULT DC CODES %

    if(dc_coeff == 0)
        huffman_dc = '00';
    elseif(dc_coeff == 1)
        huffman_dc = '010';
    elseif(dc_coeff == 2)
        huffman_dc = '011';
    elseif(dc_coeff == 3)
        huffman_dc = '100';
    elseif(dc_coeff == 4)
        huffman_dc = '101';
    elseif(dc_coeff == 5)
        huffman_dc = '110';
    elseif(dc_coeff == 6)
        huffman_dc = '1110';
    elseif(dc_coeff == 7)
        huffman_dc = '11110';
    elseif(dc_coeff == 8)
        huffman_dc = '111110';
    elseif(dc_coeff == 9)
        huffman_dc = '1111110';
    elseif(dc_coeff == 10)
        huffman_dc = '11111110';
    elseif(dc_coeff == 11)
        huffman_dc = '111111110';
    else
        huffman_dc = '00';
    end
end

function huffman_deco = Deco_DC(deco)                                      % FUNCTION TO DECODE DC COEFFICIENTS %

    if(strcmp(deco,'00'))
        huffman_deco = 0;
    elseif(strcmp(deco,'010'))
        huffman_deco = 1;
    elseif(strcmp(deco,'011'))
        huffman_deco = 2;
    elseif(strcmp(deco,'100'))
        huffman_deco = 3;
    elseif(strcmp(deco,'101'))
        huffman_deco = 4;
    elseif(strcmp(deco,'110'))
        huffman_deco = 5;
    elseif(strcmp(deco,'1110'))
        huffman_deco = 6;
    elseif(strcmp(deco,'11110'))
        huffman_deco = 7;
    elseif(strcmp(deco,'111110'))
        huffman_deco = 8;
    elseif(strcmp(deco,'1111110'))
        huffman_deco = 9;
    elseif(strcmp(deco,'11111110'))
        huffman_deco = 10;
    elseif(strcmp(deco,'111111110'))
        huffman_deco = 11;
    else
        huffman_deco = 0;
    end
end

function code = lbc_encoder(bitstream)                                      % Functions for Linear Block Codes %
    G1 = [1 1 0 1 0 0 1];
    G2 = [1 0 1 0 0 1 1];
    G3 = [1 1 1 0 1 0 0];
    c1 = and(bitstream(1),G1);
    c2 = and(bitstream(2),G2);
    c3 = and(bitstream(3),G3); 
    c4 = or(c1,c2);
    code1 = or(c3,c4);
    G11 = [1 0 0 1 1 1 0];
    G21 = [0 1 0 0 1 1 1];
    G31 = [0 0 1 1 1 0 1];
    c11 = and(bitstream(1),G11);
    c21 = and(bitstream(2),G21);
    c31 = and(bitstream(3),G31); 
    c41 = or(c11,c21);
    code = or(c31,c41);
end

function decode = lbc_decoder(cw,dbs)
    H1 = [1 0 1 1 0 0 0];
    H2 = [1 1 1 0 1 0 0];
    H3 = [1 1 0 0 0 1 0];
    H4 = [0 1 1 0 0 0 1];
    H = [1 0 1 1 0 0 0;1 1 1 0 1 0 0;1 1 0 0 0 1 0;0 1 1 0 0 0 1];
    Hd = transpose(H);
    s = cw; 
    if s == 0
        dbs = cw;
    elseif strcmp(s,'1110') == 0
        cw(1) = 0;
        dbs = cw;
    elseif strcmp(s,'0111') == 0
        cw(2) = 0;
        dbs = cw;
    elseif strcmp(s,'1101') == 0
        cw(1) = 1;
        dbs = cw;
    elseif strcmp(s,'1000') == 0
        cw(1) = 1;
        dbs = cw;
    elseif strcmp(s,'1110') == 0 
        cw(1) = 0;
        dbs = cw;
    elseif strcmp(s,'1110') == 0
        cw(1) = 0;
        dbs = cw;
    else
        cw(1) = 0;
        dbs = cw;
    end
    decode= cw + dbs;
end






