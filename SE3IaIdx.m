function w = SE3IaIdx(iS,I)
Map = [1 7 8 13 16 19; 7 2 9 14 17 20; 8 9 3 15 18 21;...
    13 14 15 4 10 11; 16 17 18 10 5 12; 19 20 21 11 12 6];
w = I(Map(iS,:));
end