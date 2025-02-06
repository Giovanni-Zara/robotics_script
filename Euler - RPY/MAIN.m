function output = controller(sequence, input, type)
    if strcmpi(type, 'euler')
        disp("EULER ANGLES");
        disp(['Sequence: ', sequence]);
        switch lower(sequence)
            case 'zyx'
                output = rotation_to_euler_zyx(input);
            case 'zxy'
                output = rotation_to_euler_zxy(input);
            case 'zyz'
                output = rotation_to_euler_zyz(input);
            case 'zxz'
                output = rotation_to_euler_zxz(input);
            case 'yxz'
                output = rotation_to_euler_yxz(input);
            case 'yxy'
                output = rotation_to_euler_yxy(input);
            case 'yzx'
                output = rotation_to_euler_yzx(input);
            case 'yzy'
                output = rotation_to_euler_yzy(input);
            case 'xyz'
                output = rotation_to_euler_xyz(input);
            case 'xzy'
                output = rotation_to_euler_xzy(input);
            case 'xzx'
                output = rotation_to_euler_xzx(input);
            case 'xyx'
                output = rotation_to_euler_xyx(input);
            otherwise
                error('Invalid sequence provided');
        end
    elseif strcmpi(type, 'rpy')
        disp("ROLL - PITCH - YAW");
        rotation_sequence = sequence(end:-1:1);
        disp(['Original Sequence: ', sequence]);
        disp(['Reversed sequence: ', rotation_sequence]);
        switch lower(rotation_sequence)
            case 'zyx'
                output = rotation_to_euler_zyx(input);
            case 'zxy'
                output = rotation_to_euler_zxy(input);
            case 'zyz'
                output = rotation_to_euler_zyz(input);
            case 'zxz'
                output = rotation_to_euler_zxz(input);
            case 'yxz'
                output = rotation_to_euler_yxz(input);
            case 'yxy'
                output = rotation_to_euler_yxy(input);
            case 'yzx'
                output = rotation_to_euler_yzx(input);
            case 'yzy'
                output = rotation_to_euler_yzy(input);
            case 'xyz'
                output = rotation_to_euler_xyz(input);
            case 'xzy'
                output = rotation_to_euler_xzy(input);
            case 'xzx'
                output = rotation_to_euler_xzx(input);
            case 'xyx'
                output = rotation_to_euler_xyx(input);
            otherwise
                error('Invalid sequence provided');
        end
        %disp(["output from function", output]);
        disp("IMPORTANTE -- CASO RPY DEVI SOSTITUIRE PRIMO ANGOLO CON TERZO");
        disp("!!! sostituzione va eseguita sia per caso positivo che per caso negativo !!!")
        disp("!! riportata qui sotto UNA delle due già invertita, l'altra è da invertire a mano !!")
        output([1, 3]) = output([3, 1]);
        disp(["output AFTER SWTICH", output]);

    else
        error('Invalid type of angle sequence provided');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


R = [sqrt(6)/4   sqrt(2)/4   -sqrt(2)/2;
     -sqrt(6)/4  -sqrt(2)/4  -sqrt(2)/2;
     -1/2        sqrt(3)/2    0];

a = controller('zyx', R, 'euler');