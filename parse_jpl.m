function [E_Sun,VE_Sun,E_Earth,VE_Earth,E_Jupiter,VE_Jupiter,E_Bennu,VE_Bennu,E_OR,VE_OR] = parse_jpl()
%Parses JPL ephemeris files and places values in arrays
addpath('Ephemerides');

[E_Sun,VE_Sun] = parse_file('Sun.txt');
[E_Earth,VE_Earth] = parse_file('Earth.rtf');
[E_Jupiter,VE_Jupiter] = parse_file('Jupiter.rtf');
[E_Bennu,VE_Bennu] = parse_file('Bennu.rtf');
[E_OR,VE_OR] = parse_file('OR.rtf');

end

function [r,v] = parse_file(file)
%Given file name, parses file and returns array
fileID = fopen(file);
foo = textscan(fileID,'%s',1,'delimiter','\n');
k = 1;
%Goes through file until ephemeris data starts
while strcmp(foo{1}{1},'$$SOE\') ~= 1
    foo = textscan(fileID,'%s',1,'delimiter','\n');
    k = k + 1;
end
r = [];
v = [];
%Continues through file until end of ephemeris data
line = textscan(fileID,'%s',4,'delimiter','\n');
while strcmp(line{1}{1},'$$EOE\') ~= 1
    %Making sure there is a space after each = sign
    for k = 2:3
        bas = strfind(line{1,1}{k},'=');
        bar = strfind(line{1,1}{k},'-');
        for j = 1:length(bas)
            for i = 1:length(bar)
                if bar(i) == bas(j) + 1
                    foo = insertBefore(line{1,1}{k},bar(i),' ');
                    line{1,1}{k} = foo;
                    %fprintf('i is %i and j is %i.\n',bar(i),bas(j));
                    %display(line{1,1}{k});
                    %When a space is added, update the new positions of
                    %later = and - 
                    for l = 1:length(bas)
                        if bas(l) > bas(j)
                            bas(l) = bas(l) + 1;
                        end
                    end
                    for l = 1:length(bar)
                        if bar(l) > bar(i)
                            bar(l) = bar(l) + 1;
                        end
                    end
                end
            end
        end
    end
    
    r_line = strsplit(line{1}{2});
    v_line = strsplit(line{1}{3});
    %display(r_line);
    %display(v_line);
    
    rl = strlength(r_line{1,9});
    vl = strlength(v_line{1,6});
    
    r = [r; str2double(r_line{1,3}), str2double(r_line{1,6}), str2double(extractBefore(r_line{1,9},rl))]; %km
    v = [v; str2double(v_line{1,2}), str2double(v_line{1,4}), str2double(extractBefore(v_line{1,6},vl))]; %km/s
    
    line = textscan(fileID,'%s',4,'delimiter','\n');
    
    k = k + 1;
end
fprintf('%s JPL ephem file parse finished.\n',file);
end