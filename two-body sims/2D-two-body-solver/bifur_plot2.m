function bifur_plot2(csv_directory, lambdaStop, logLambdaAxis)
% lambdaStop is the value at which you want to stop plotting coupling.
% logLambdaAxis is a boolean that, if true, will plot coupling on a log axis.

    % Initialize the figure
    lambda_plot = figure('Name', 'Lambda-visualized Bifurcation Diagram');
    hold on
    grid on

    % Make the graph look pretty and informative
    set(gca, 'Layer', 'Top');
    set(gcf, 'color', 'w');
    title(csv_directory(end-3:end));
    xlabel('u_1(0)')
    ylabel('u_2(0)')
    zlabel('\lambda')
    view(3);
    
    if logLambdaAxis
        set(gca,'ZScale','log');
    end

    % Parse csv's and plot paths
    file_attribs = dir(csv_directory);
    for i = 1:length(file_attribs)
        if ~file_attribs(i).isdir && file_attribs(i).name(1:9) ~= "norender_"
            bifur_path = path_from_csv([csv_directory,'/',file_attribs(i).name], lambdaStop);
            draw_path(bifur_path);
        end
    end

    hold off
    
end

function output = draw_path(bifur_path)
    %{ 
    This function generates a 3D plot of the given paths using patch().
    Input:
        bifur_path  - a path structs containing the row vectors lambda, p0(0), p1(0), p2(0),
                      and some other stuff like indices of the changes in conjugate
                      points, links to reference images, etc.
    Output:
        output      - a patch handle list
    %}

    % Draw a patch with a different line style each time the number of conjugate points
    % changes
    osi = bifur_path.osi;
    output = zeros(1,length(osi));
    for i = 1:(length(osi)-1)
        % Define line style based on number of conjugate points
        if bifur_path.num_conj_pt(osi(i)) >= 2
            color = 'black';
        elseif bifur_path.num_conj_pt(osi(i)) == 1
            color = 'red';
        else
            color = 'blue';
        end
        
        p3 = bifur_path.p3(osi(i):osi(i+1));
        p6 = bifur_path.p6(osi(i):osi(i+1));
        lambda = bifur_path.lambda(osi(i):osi(i+1));
        u1 = ((1+lambda).*p3+lambda.*p6)./(2*lambda+1);
        u2 = (lambda.*p3+(1+lambda).*p6)./(2*lambda+1);
        if osi(i+1)-osi(i) == 1
            pat = patch(u1, u2, lambda, 'black', 'EdgeColor', color, ...
                'LineWidth', 1);
        else
            pat = patch([u1 nan], [u2 nan], [lambda nan], 'black', ...
                'EdgeColor', color, 'LineWidth', 1);
        end
        
        % Set some properties of the patch which we'll use later
        pat.DisplayName = bifur_path.name;
        set(pat, 'EdgeLighting', 'none');        
        
        % Give the patch user data which is confusing and bad and oh god no please
        pat.UserData = bifur_path.num_conj_pt(osi(i));
        
        output(i) = pat;
    end

end

function path_struct = path_from_csv(filename, lambdaStop)
    %{
    Creates a path struct from the given csv file with the format:
            lambda,p1(0),p2(0),p3(0),p4(0),p5(0),p6(0),num_conj_pt
    and a header of the format:
            j                                           <--- number of reference images
            image_path_1,p1_1,p2_1,p3_1,p4_1,p5_1,p6_1  -+    
            image_path_2,p1_2,p2_2,p3_2,p4_2,p5_2,p6_2   |   path to reference image followed
            ...                                          |    by the location to plot at
            image_path_j,p1_j,p2_j,p3_j,p4_j,p5_j,p6_j  -+          
    Input:
        filename    - path to csv file to parse.   
    Output:
        path_struct - a struct representing a path in the
                        bifurcation plot
    %}
    
    % Get header entries
    fd = fopen(filename, 'r');
    num_entries = str2double(fgetl(fd));
    fclose(fd);

    % Read csv entries into a matrix
    A = readmatrix(filename, 'NumHeaderLines', num_entries+1);
    A = sortrows(A);

    % Create struct fields
    path_struct.name = replace(filename,"_","\_");
    
    stopIdx = find(A(:,1)==lambdaStop);
    if isempty(stopIdx)
        stopIdx = length(A(:,1));
    end
    path_struct.lambda = A(1:stopIdx,1)';
    path_struct.p1 = A(1:stopIdx,2)';
    path_struct.p2 = A(1:stopIdx,3)';
    path_struct.p3 = A(1:stopIdx,4)'; 
    path_struct.p4 = A(1:stopIdx,5)'; 
    path_struct.p5 = A(1:stopIdx,6)'; 
    path_struct.p6 = A(1:stopIdx,7)'; 
    path_struct.num_conj_pt = A(1:stopIdx,8)';
    
    % Find OSI
    num_conj_diff = path_struct.num_conj_pt(2:end) - path_struct.num_conj_pt(1:(end-1));
    path_struct.osi(1) = 1; % Ensure start point is 1
    for i = 1:length(num_conj_diff)
        if num_conj_diff(i) ~= 0
            path_struct.osi(end+1) = i+1;
        end
    end
    path_struct.osi(end+1) = length(path_struct.lambda); % Ensure endpoint

end