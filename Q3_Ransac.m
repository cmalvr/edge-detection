%Assignment 1 - COMP558
%QUESTION 3 - RANSAC
%Author: Camila Alvarez
%Date: 2021 October 3rd

makeEdgeList.m; %running script to return edges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% A)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges1 = edges;
M = 6000;
distance = 10; %threshold distance
anglethreshold = 0.1; %threshold angle
cmin = 100; %minimum c
line_model = [];

   for i = 1:M
    %Count c lines that are within distance and angle threshold (i.e) on the consensus set           
       
     %drawing datasample line
       rnum = randi (size(edges1,1));
       s_edge = edges1(rnum,:);
       m = s_edge(4)/s_edge(3); %slope = tangent
       y = (m * (0-s_edge(1))) + edge(2); %point slope form
       consensus = []; %initializing matriz to store points belonging to consensus set
       deletion = []; %initializing matriz to store
       % ITERATING THROUGH N-1 ELEMENTS
   
       for j = 1:(size(edges1,1))
            if(j == rnum)
                continue
            end      
            t_edge= edges1(j,:);

          %obtaining angle difference
            angle = abs(acos(s_edge(3))-acos(t_edge(3)));
    
            %obtaining distance from point in t_edge(1),t_edge(2) to line defined by (0, y) and (s_edge(1),s_edge(2))
            d = pointdistance(t_edge(1),t_edge(2),0,y,s_edge(1),s_edge(2));
          %Now we check if edge its within thresholds
            
            if (angle <= anglethreshold) && (d <= distance)
              consensus = [consensus; t_edge(1) t_edge(2) acos(t_edge(3))]; %add row to array x y angle index in edge   
              deletion = [deletion; t_edge];
            end   
       end %FINISH ITERATING THROUGH N-1 ELEMENTS
        %If size of consensus set is larger than cmin
       if (length(consensus) >= cmin)  
           consensus = [consensus;s_edge(1) s_edge(2) acos(s_edge(3))]; 
           deletion = [deletion; s_edge];
             [row, col] =  size(consensus);
             [row_d, col_d] =  size(deletion);
           %Iterating through consensus set to get average values
           av_x = 0;
           av_y = 0;
           av_a = 0;

           for h=1:row        
               av_x = av_x + consensus(h,1);
               av_y = av_y + consensus(h,2);
               av_a = av_a + consensus(h,3); 
           end

           for t = 1:row_d
              
              C1 = edges1(:,1) == deletion(t,1);
              C2 = edges1(:,2) == deletion(t,2);
              C3 = edges1(:,3) == deletion(t,3);
              C4 = edges1(:,4) == deletion(t,4);
              F = C1 & C2 & C3 & C4;
              edges1 (F,:) = [];
           end

           lol = edges1;
            %Computing average values to construct line model
           av_x = av_x/length(consensus);
           av_y = av_y/length(consensus);
           av_a = av_a/length(consensus);
           line_model = [line_model; av_x av_y av_a];
           edges1 = [edges1; av_x av_y cos(av_a) sin(av_y)];
       end
    
   end % END OF RANSAC

figure('Name','Q3) Edge Detection');
imshow(im)
hold on

start = 1;
[row_1, col_1] = size(line_model);
if (row_1)>20
    start = row_1 -10;
end

%Plotting 10 last line models since those are expected to be more accurate 
for k=start:row_1
    x0 = line_model(k,1);
    y0 = line_model(k,2);
    angle = line_model(k,3);
    x = 0 : 250;
    y = tan(angle)*(x-x0) + y0;
    plot(x,y,'r')
end
%Function used to compute point to line passing through two points distance
function d = pointdistance(x0,y0,x1,y1,x2,y2)

  d = abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1))./ sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);


end

