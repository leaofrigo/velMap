function letter = num2letter(i)
if and(rem(i,1)==0,i>0) ==0
    error('Number must be integer bigger then zero smaller then 18278')
end
if i>18278
    error('Number must be integer bigger then zero smaller then 18278')
end
j=0;
k=0;
letter2 = '';
letter3 = '';
while i>26
    i = i-26;
    j=j+1;
end
while j > 26
    j = j-26;
    k=k+1;
end
if i == 1
    letter1 = 'A';
elseif i==2
    letter1 = 'B';
elseif i==3
    letter1 = 'C';
elseif i==4
    letter1 = 'D';
elseif i==5
    letter1 = 'E';
elseif i==6
    letter1 = 'F';
elseif i==7
    letter1 = 'G';
elseif i==8
    letter1 = 'H';
elseif i==9
    letter1 = 'I';
elseif i==10
    letter1 = 'J';
elseif i==11
    letter1 = 'K';
elseif i==12
    letter1 = 'L';
elseif i==13
    letter1 = 'M';
elseif i==14
    letter1 = 'N';
elseif i==15
    letter1 = 'O';
elseif i==16
    letter1 = 'P';
elseif i==17
    letter1 = 'Q';
elseif i==18
    letter1 = 'R';
elseif i==19
    letter1 = 'S';
elseif i==20
    letter1 = 'T';
elseif i==21
    letter1 = 'U';
elseif i==22
    letter1 = 'V';
elseif i==23
    letter1 = 'W';
elseif i==24
    letter1 = 'X';
elseif i==25
    letter1 = 'Y';
elseif i==26
    letter1 = 'Z';
end
if j == 1
    letter2 = 'A';
elseif j == 2
    letter2 = 'B';
elseif j == 3
    letter2 = 'C';
elseif j == 4
    letter2 = 'D';
elseif j == 5
    letter2 = 'E';
elseif j == 6
    letter2 = 'F';
elseif j == 7
    letter2 = 'G';
elseif j == 8
    letter2 = 'H';
elseif j == 9
    letter2 = 'I';
elseif j == 10
    letter2 = 'J';
elseif j == 11
    letter2 = 'K';
elseif j == 12
    letter2 = 'L';
elseif j == 13
    letter2 = 'M';
elseif j == 14
    letter2 = 'N';
elseif j == 15
    letter2 = 'O';
elseif j == 16
    letter2 = 'P';
elseif j == 17
    letter2 = 'Q';
elseif j == 18
    letter2 = 'R';
elseif j == 19
    letter2 = 'S';
elseif j == 20
    letter2 = 'T';
elseif j == 21
    letter2 = 'U';
elseif j == 22
    letter2 = 'V';
elseif j == 23
    letter2 = 'W';
elseif j == 24
    letter2 = 'X';
elseif j == 25
    letter2 = 'Y';
elseif j == 26
    letter2 = 'Z';
end
if k == 1
    letter3 = 'A';
elseif k == 2
    letter3 = 'B';
elseif k == 3
    letter3 = 'C';
elseif k == 4
    letter3 = 'D';
elseif k == 5
    letter3 = 'E';
elseif k == 6
    letter3 = 'F';
elseif k == 7
    letter3 = 'G';
elseif k == 8
    letter3 = 'H';
elseif k == 9
    letter3 = 'I';
elseif k == 10
    letter3 = 'J';
elseif k == 11
    letter3 = 'K';
elseif k == 12
    letter3 = 'L';
elseif k == 13
    letter3 = 'M';
elseif k == 14
    letter3 = 'N';
elseif k == 15
    letter3 = 'O';
elseif k == 16
    letter3 = 'P';
elseif k == 17
    letter3 = 'Q';
elseif k == 18
    letter3 = 'R';
elseif k == 19
    letter3 = 'S';
elseif k == 20
    letter3 = 'T';
elseif k == 21
    letter3 = 'U';
elseif k == 22
    letter3 = 'V';
elseif k == 23
    letter3 = 'W';
elseif k == 24
    letter3 = 'X';
elseif k == 25
    letter3 = 'Y';
elseif k == 26
    letter3 = 'Z';
end
letter = [letter3 letter2 letter1];