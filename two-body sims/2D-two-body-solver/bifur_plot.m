function bifur_plot(csv_directory)

    % Initialize the figure
    fig = figure('Name', 'Bifurcation Diagram (p_1,p_2,p_3)');
    hold on
    grid on
    %axis equal
   
    % Make the graph look pretty and informative
    xlabel('p_1(0)')
    ylabel('p_2(0)')
    zlabel('p_3(0)')

    colormap(copper);
    c = colorbar;
    c.Label.String = 'mu';
    daspect([1 1 1]);
    set(gca, 'Layer', 'Top');
    view(3);

    % Parse csv's and plot paths
    file_attribs = dir(csv_directory);
    bifur_paths= [];
    for i = 1:length(file_attribs)
        if ~file_attribs(i).isdir && file_attribs(i).name(1:9) ~= "norender_"
            bifur_path = path_from_csv([csv_directory,'/',file_attribs(i).name]);
            draw_path(bifur_path, 'default');
            bifur_paths = [bifur_paths bifur_path];
        end
    end

    % Set a custom datatip for the paths
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@(~,event) update_data_tip(event));
    
    hold off
    
    % Initialize the p4+ figure
    fig = figure('Name', 'Bifurcation Diagram (p_4,p_5,p_6)');
    hold on
    grid on
    %axis equal
    
    % Make the graph look pretty and informative
    xlabel('p_4(0)')
    ylabel('p_5(0)')
    zlabel('p_6(0)')

    colormap(copper);
    c = colorbar;
    c.Label.String = 'mu';
    daspect([1 1 1]);
    set(gca, 'Layer', 'Top');
    view(3);

    % Parse csv's and plot paths
    file_attribs = dir(csv_directory);
    for i = 1:length(file_attribs)
        if ~file_attribs(i).isdir && file_attribs(i).name(1:9) ~= "norender_"
            bifur_path = path_from_csv([csv_directory,'/',file_attribs(i).name]);
            draw_path(bifur_path, 'p4+');
        end
    end

    % Set a custom datatip for the paths
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@(~,event) update_data_tip(event));
    
    hold off
    

    % Initialize the second figure
    mu_plot = figure('Name', 'Mu-visualized Bifurcation Diagram');
    hold on
    grid on

    % Make the graph look pretty and informative
    set(gca, 'Layer', 'Top');
    xlabel('p_1(0)')
    ylabel('p_2(0)')
    zlabel('mu')
    colormap(copper);
    c = colorbar;
    c.Label.String = 'p3';
    daspect([1 1 1/100]);
    view(3);

    % Parse csv's and plot paths
    file_attribs = dir(csv_directory);
    for i = 1:length(file_attribs)
        if ~file_attribs(i).isdir && file_attribs(i).name(1:9) ~= "norender_"
            bifur_path = path_from_csv([csv_directory,'/',file_attribs(i).name]);
            draw_path(bifur_path, 'mu');
        end
    end

    hold off

    % bifur_graph(bifur_paths);

end

function output = draw_path(bifur_path, tag)
    %{ 
    This function generates a 3D plot of the given paths using patch().
    Input:
        bifur_path  - a path structs containing the row vectors mu, p0(0), p1(0), p2(0),
                      and some other stuff like indices of the changes in conjugate
                      points, links to reference images, etc.
        tag         - char vector with a tag saying to flip axes or not
    Output:
        output      - a patch handle list
    %}

    kPathClickCallback = 1;
    kPathStartDisabled = 1;

    % Draw a patch with a different line style each time the number of conjugate points
    % changes
    osi = bifur_path.osi;
    output = zeros(1,length(osi));
    imgIdx= 1;
    for i = 1:(length(osi)-1)
        % Define line style based on number of conjugate points
        if bifur_path.num_conj_pt(osi(i)) >= 2
            line_style = ':';
        elseif bifur_path.num_conj_pt(osi(i)) == 1
            line_style = '--';
        else
            line_style = '-';
        end
        
        % Render based on what type of graph we're making
        switch tag
        case 'default'
            if osi(i+1)-osi(i) == 1
                pat = patch(bifur_path.p_mat(osi(i):osi(i+1),1)', ...
                    bifur_path.p_mat(osi(i):osi(i+1),2)', ...
                    bifur_path.p_mat(osi(i):osi(i+1),3)', ...
                    bifur_path.mu(osi(i):osi(i+1)), ...
                    'LineStyle', line_style,'FaceColor', 'none', 'EdgeColor', 'interp');
            else
                pat = patch([bifur_path.p_mat(osi(i):osi(i+1),1)' nan], ...
                    [bifur_path.p_mat(osi(i):osi(i+1),2)' nan], ...
                    [bifur_path.p_mat(osi(i):osi(i+1),3)' nan], ...
                    [bifur_path.mu(osi(i):osi(i+1)) nan], ...
                    'LineStyle', line_style,'FaceColor', 'none', 'EdgeColor', 'interp');
            end
        case 'mu'
            if osi(i+1)-osi(i) == 1
                pat = patch(bifur_path.p_mat(osi(i):osi(i+1),1)', ...
                    bifur_path.p_mat(osi(i):osi(i+1),2)', ...
                    bifur_path.mu(osi(i):osi(i+1)), ...
                    bifur_path.p_mat(osi(i):osi(i+1),3)', ...
                    'LineStyle', line_style,'FaceColor', 'none', 'EdgeColor', 'interp');
            else
                pat = patch([bifur_path.p_mat(osi(i):osi(i+1),1)' nan], ...
                    [bifur_path.p_mat(osi(i):osi(i+1),2)' nan], ...
                    [bifur_path.mu(osi(i):osi(i+1)) nan], ...
                    [bifur_path.p_mat(osi(i):osi(i+1),3)' nan], ...
                    'LineStyle', line_style,'FaceColor', 'none', 'EdgeColor', 'interp');
            end
        case 'p4+'
            if osi(i+1)-osi(i) == 1
                pat = patch(bifur_path.p_mat(osi(i):osi(i+1),4)', ...
                    bifur_path.p_mat(osi(i):osi(i+1),5)', ...
                    bifur_path.p_mat(osi(i):osi(i+1),6)', ...
                    bifur_path.mu(osi(i):osi(i+1)), ...
                    'LineStyle', line_style,'FaceColor', 'none', 'EdgeColor', 'interp');
            else
                pat = patch([bifur_path.p_mat(osi(i):osi(i+1),4)' nan], ...
                    [bifur_path.p_mat(osi(i):osi(i+1),5)' nan], ...
                    [bifur_path.p_mat(osi(i):osi(i+1),6)' nan], ...
                    [bifur_path.mu(osi(i):osi(i+1)) nan], ...
                    'LineStyle', line_style,'FaceColor', 'none', 'EdgeColor', 'interp');
            end
        end
        
        % Draw reference image for this segment
        if i==1 || i==length(osi)-1 || i==ceil((length(osi)-1)/2)
            image_handle = draw_image_surface(bifur_path.image_paths{imgIdx}, ...
                bifur_path.image_pts(imgIdx,:), tag);
            imgIdx= imgIdx+1;
        end
        
        % Set some properties of the patch which we'll use later
        pat.DisplayName = bifur_path.name;
        set(pat, 'EdgeLighting', 'none');
        if kPathClickCallback
            set(pat, 'ButtonDownFcn', @on_path_click);
        end
        if kPathStartDisabled
            pat.EdgeAlpha = 0.2;
            if image_handle ~= -1
                set(image_handle, 'Visible', 'off');
            end
        end        
        
        % Give the patch user data which is confusing and bad and oh god no please
        pat.UserData = {image_handle; bifur_path.num_conj_pt(osi(i))};
        
        output(i) = pat;
    end

    % Draw a dot at mu = 0 if the path contains mu=0
    if bifur_path.mu(1) < 0.001
        switch tag
        case 'default'
            plot3(bifur_path.p_mat(1,1), bifur_path.p_mat(1,2), bifur_path.p_mat(1,3), 'ko');
        case 'mu'
            plot3(bifur_path.p_mat(1,1), bifur_path.p_mat(1,2), bifur_path.mu(1), 'ko');
        case 'p4+'
            plot3(bifur_path.p_mat(1,4), bifur_path.p_mat(1,5), bifur_path.p_mat(1,6), 'ko');
        end
    end
end

function path_struct = path_from_csv(filename)
    %{
    Creates a path struct from the given csv file with the format:
            mu,p1(0),p2(0),p3(0),p4(0),p5(0),p6(0),num_conj_pt
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
    num_entries = str2num(fgetl(fd));

    % Make a list of image paths and render locations
    image_paths = cell(num_entries,1);
    image_pts = zeros(num_entries,3);
    for i = 1:num_entries
        split_line = split(fgetl(fd), ',');
        image_paths(i) = split_line(1);
        image_pts(i,:) = str2num(char(split_line(2:4)));
    end
    fclose(fd);

    % Read csv entries into a matrix
    A = readmatrix(filename, 'NumHeaderLines', num_entries+1);
    A = sortrows(A);

    % Create struct fields
    path_struct.name = replace(filename,"_","\_");
    path_struct.image_paths = image_paths;
    path_struct.image_pts = image_pts;
    path_struct.p_mat = A(:,2:7);
    path_struct.mu = A(:,1)';
    path_struct.p1 = A(:,2)';
    path_struct.p2 = A(:,3)';
    path_struct.p3 = A(:,4)'; 
    path_struct.p4 = A(:,5)'; 
    path_struct.p5 = A(:,6)'; 
    path_struct.p6 = A(:,7)'; 
    path_struct.num_conj_pt = A(:,8)';

    % Find OSI
    num_conj_diff = path_struct.num_conj_pt(2:end) - path_struct.num_conj_pt(1:(end-1));
    path_struct.osi(1) = 1; % Ensure start point is 1
    for i = 1:length(num_conj_diff)
        if num_conj_diff(i) ~= 0
            path_struct.osi(end+1) = i+1;
        end
    end
    path_struct.osi(end+1) = length(path_struct.mu); % Ensure endpoint

end

function output_text = update_data_tip(event_obj)
    %{
    Callback function for MatLab's data tip. Displays a tooltip with a point's
    number of conjugate points, p_1 through p_3, and mu values.
    Input:
        event_obj - an event struct containing the handle of the object selected,
                    the position of the mouse, and the index of the data point
    %}
    target = event_obj.Target;
    position = event_obj.Position;
    index = event_obj.DataIndex;

    if target.DisplayName == ""
        output_text = { 
                        ['num conj pts: ', num2str(target.UserData{2})],...
                        ['p_1(0): ', num2str(position(1))],...
                        ['p_2(0): ', num2str(position(2))],...
                        ['p_3(0): ', num2str(position(3))],...
                        ['mu: ', num2str(target.CData(index))] };
    elseif target.DisplayName(1:6) == "image_"
        output_text = { target.DisplayName }; 
    else
        output_text = { target.DisplayName,...
                        ['num conj pts: ', num2str(target.UserData{2})],...
                        ['p_1(0): ', num2str(position(1))],...
                        ['p_2(0): ', num2str(position(2))],...
                        ['p_3(0): ', num2str(position(3))],...
                        ['mu: ', num2str(target.CData(index))] };
    end

end

function on_image_click(surf_handle, ~) 
    %{
    Callback function executed when an image is clicked in the figure
    which fades the image considerably.
    NOTE: this callback is only executed if kImageClickCallback is true in
    draw_image_surface, below
    Input:
        surf_handle - the handle of the surface object clicked
    %} 

    other_alphas = surf_handle.UserData;
    surf_handle.UserData = surf_handle.AlphaData;
    surf_handle.AlphaData = other_alphas;
end

function on_path_click(patch_handle, ~)
    image_handle = patch_handle.UserData{1};
    if image_handle == -1
        return;
    end

    if patch_handle.EdgeAlpha < 1
        patch_handle.EdgeAlpha = 1;
        set(image_handle, 'Visible', 'on');
    else
        patch_handle.EdgeAlpha = 0.2;
        set(image_handle, 'Visible', 'off');
    end
end

function output = draw_image_surface(image_path, image_point, tag)
    %{
    Draws an image near the given point as a 3d surface
    Input:
        image_path  - path to image to render
        image_point - a vector [p1 p2 p3] which is the point near which
                      the image will be drawn 
        tag         - a character vector indicating whether to switch axes or not
    %}

    kImageClickCallback = 0;
    kImageDim = 10;
    kImageZOffset = 7;

    % Get image data
    [I,map,transparency] = imread(image_path);

    % Create surface
    range_x = (image_point(1) - kImageDim):kImageDim:(image_point(1) + kImageDim); % I made the step size the image dimension
    range_y = (image_point(2) - kImageDim):kImageDim:(image_point(2) + kImageDim); %   so that the polygons would be at a min
    [X, Y] = meshgrid(range_x, range_y);
    switch tag
    case 'default'
        Z = (image_point(3) + kImageZOffset) * ones(length(X(:,1)), length(X(1,:)));
    case 'mu'
        Z = (0) * ones(length(X(:,1)), length(X(1,:)));
    case 'p4+'
        output = -1;
        return
    end
    
    % Render
    tex = warp(X, Y, Z, I);
    
    % Set properties of the surface
    tex.DisplayName = strcat("image_",replace(image_path,"_","\_"));
    set(tex, 'AlphaData', transparency);
    set(tex, 'FaceAlpha', 'texturemap');

    % If applicable, set the click callback function
    if kImageClickCallback
        set(tex, 'UserData', transparency * 0.1); % Store a copy of the transparency matrix in UserData
        set(tex, 'ButtonDownFcn', @on_image_click); % Define click callback for the image surface
    end
    set(tex, 'FaceLighting', 'none');

    output = tex;
end
