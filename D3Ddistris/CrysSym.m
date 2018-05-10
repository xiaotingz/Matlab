function [O] = CrysSym

O = zeros(3,3,24);
    
% !     1       h  k  l
	O(1,1,1)=1.0;
	O(1,2,1)=0.0;
	O(1,3,1)=0.0;
	O(2,1,1)=0.0;
	O(2,2,1)=1.0;
	O(2,3,1)=0.0;
	O(3,1,1)=0.0;
	O(3,2,1)=0.0;
	O(3,3,1)=1.0;
% c  2		 l  h  k
	O(1,1,2)=0.0;
	O(1,2,2)=0.0;
	O(1,3,2)=1.0;
	O(2,1,2)=1.0;
	O(2,2,2)=0.0;
	O(2,3,2)=0.0;
	O(3,1,2)=0.0;
	O(3,2,2)=1.0;
	O(3,3,2)=0.0;
% c  3          k  l  h
	O(1,1,3)=0.0;
	O(1,2,3)=1.0;
	O(1,3,3)=0.0;
	O(2,1,3)=0.0;
	O(2,2,3)=0.0;
	O(2,3,3)=1.0;
	O(3,1,3)=1.0;
	O(3,2,3)=0.0;
	O(3,3,3)=0.0;
% c  4	  -k l  -h
	O(1,1,4)=0.0;
	O(1,2,4)=-1.0;
	O(1,3,4)=0.0;
	O(2,1,4)=0.0;
	O(2,2,4)=0.0;
	O(2,3,4)=1.0;
	O(3,1,4)=-1.0;
	O(3,2,4)=0.0;
	O(3,3,4)=0.0;
% c  5       -k -l h 
	O(1,1,5)=0.0;
	O(1,2,5)=-1.0;
	O(1,3,5)=0.0;
	O(2,1,5)=0.0;
	O(2,2,5)=0.0;
	O(2,3,5)=-1.0;
	O(3,1,5)=1.0;
	O(3,2,5)=0.0;
	O(3,3,5)=0.0;
% c  6	  k  -l -h
	O(1,1,6)=0.0;
	O(1,2,6)=1.0;
	O(1,3,6)=0.0;
	O(2,1,6)=0.0;
	O(2,2,6)=0.0;
	O(2,3,6)=-1.0;
	O(3,1,6)=-1.0;
	O(3,2,6)=0.0;
	O(3,3,6)=0.0;
% c  7     -l h  -k 
	O(1,1,7)=0.0;
	O(1,2,7)=0.0;
	O(1,3,7)=-1.0;
	O(2,1,7)=1.0;
	O(2,2,7)=0.0;
	O(2,3,7)=0.0;
	O(3,1,7)=0.0;
	O(3,2,7)=-1.0;
	O(3,3,7)=0.0;
% c  8      -l -h k
	O(1,1,8)=0.0;
	O(1,2,8)=0.0;
	O(1,3,8)=-1.0;
	O(2,1,8)=-1.0;
	O(2,2,8)=0.0;
	O(2,3,8)=0.0;
	O(3,1,8)=0.0;
	O(3,2,8)=1.0;
	O(3,3,8)=0.0;
% c    9       l  -h -k
	O(1,1,9)=0.0;
	O(1,2,9)=0.0;
	O(1,3,9)=1.0;
	O(2,1,9)=-1.0;
	O(2,2,9)=0.0;
	O(2,3,9)=0.0;
	O(3,1,9)=0.0;
	O(3,2,9)=-1.0;
	O(3,3,9)=0.0;
% c 10		 -h k  -l 
	O(1,1,10)=-1.0;
	O(1,2,10)=0.0;
	O(1,3,10)=0.0;
	O(2,1,10)=0.0;
	O(2,2,10)=1.0;
	O(2,3,10)=0.0;
	O(3,1,10)=0.0;
	O(3,2,10)=0.0;
	O(3,3,10)=-1.0;
% c 11         -h -k l
	O(1,1,11)=-1.0;
	O(1,2,11)=0.0;
	O(1,3,11)=0.0;
	O(2,1,11)=0.0;
	O(2,2,11)=-1.0;
	O(2,3,11)=0.0;
	O(3,1,11)=0.0;
	O(3,2,11)=0.0;
	O(3,3,11)=1.0;
% c 12	      h  -k -l 
	O(1,1,12)=1.0;
	O(1,2,12)=0.0;
	O(1,3,12)=0.0;
	O(2,1,12)=0.0;
	O(2,2,12)=-1.0;
	O(2,3,12)=0.0;
	O(3,1,12)=0.0;
	O(3,2,12)=0.0;
	O(3,3,12)=-1.0;
% c 13           -l -k -h 
	O(1,1,13)=0.0;
	O(1,2,13)=0.0;
	O(1,3,13)=-1.0;
	O(2,1,13)=0.0;
	O(2,2,13)=-1.0;
	O(2,3,13)=0.0;
	O(3,1,13)=-1.0;
	O(3,2,13)=0.0;
	O(3,3,13)=0.0;
% c 14	       l  -k h  
	O(1,1,14)=0.0;
	O(1,2,14)=0.0;
	O(1,3,14)=1.0;
	O(2,1,14)=0.0;
	O(2,2,14)=-1.0;
	O(2,3,14)=0.0;
	O(3,1,14)=1.0;
	O(3,2,14)=0.0;
	O(3,3,14)=0.0;
% c 15           l  k  -h 
	O(1,1,15)=0.0;
	O(1,2,15)=0.0;
	O(1,3,15)=1.0;
	O(2,1,15)=0.0;
	O(2,2,15)=1.0;
	O(2,3,15)=0.0;
	O(3,1,15)=-1.0;
	O(3,2,15)=0.0;
	O(3,3,15)=0.0;
% c 16           -l k  h 
	O(1,1,16)=0.0;
	O(1,2,16)=0.0;
	O(1,3,16)=-1.0;
	O(2,1,16)=0.0;
	O(2,2,16)=1.0;
	O(2,3,16)=0.0;
	O(3,1,16)=1.0;
	O(3,2,16)=0.0;
	O(3,3,16)=0.0;
% c 17            -h -l -k 
	O(1,1,17)=-1.0;
	O(1,2,17)=0.0;
	O(1,3,17)=0.0;
	O(2,1,17)=0.0;
	O(2,2,17)=0.0;
	O(2,3,17)=-1.0;
	O(3,1,17)=0.0;
	O(3,2,17)=-1.0;
	O(3,3,17)=0.0;
% c 18		    h  -l k  
	O(1,1,18)=1.0;
	O(1,2,18)=0.0;
	O(1,3,18)=0.0;
	O(2,1,18)=0.0;
	O(2,2,18)=0.0;
	O(2,3,18)=-1.0;
	O(3,1,18)=0.0;
	O(3,2,18)=1.0;
	O(3,3,18)=0.0;
% c 19          h  l  -k
	O(1,1,19)=1.0;
	O(1,2,19)=0.0;
	O(1,3,19)=0.0;
	O(2,1,19)=0.0;
	O(2,2,19)=0.0;
	O(2,3,19)=1.0;
	O(3,1,19)=0.0;
	O(3,2,19)=-1.0;
	O(3,3,19)=0.0;
% c 20	       -h l  k
	O(1,1,20)=-1.0;
	O(1,2,20)=0.0;
	O(1,3,20)=0.0;
	O(2,1,20)=0.0;
	O(2,2,20)=0.0;
	O(2,3,20)=1.0;
	O(3,1,20)=0.0;
	O(3,2,20)=1.0;
	O(3,3,20)=0.0;
% c 21          -k -h -l 
	O(1,1,21)=0.0;
	O(1,2,21)=-1.0;
	O(1,3,21)=0.0;
	O(2,1,21)=-1.0;
	O(2,2,21)=0.0;
	O(2,3,21)=0.0;
	O(3,1,21)=0.0;
	O(3,2,21)=0.0;
	O(3,3,21)=-1.0;
% c 22	      k  -h l
	O(1,1,22)=0.0;
	O(1,2,22)=1.0;
	O(1,3,22)=0.0;
	O(2,1,22)=-1.0;
	O(2,2,22)=0.0;
	O(2,3,22)=0.0;
	O(3,1,22)=0.0;
	O(3,2,22)=0.0;
	O(3,3,22)=1.0;  
% c 23          k  h  -l
	O(1,1,23)=0.0;
	O(1,2,23)=1.0;
	O(1,3,23)=0.0;
	O(2,1,23)=1.0;
	O(2,2,23)=0.0;
	O(2,3,23)=0.0;
	O(3,1,23)=0.0;
	O(3,2,23)=0.0;
	O(3,3,23)=-1.0;
% c 24           -k h  l
	O(1,1,24)=0.0;
	O(1,2,24)=-1.0;
	O(1,3,24)=0.0;
	O(2,1,24)=1.0;
	O(2,2,24)=0.0;
	O(2,3,24)=0.0;
	O(3,1,24)=0.0;
	O(3,2,24)=0.0;
	O(3,3,24)=1.0;
    
end