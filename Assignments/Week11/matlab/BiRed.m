function [ A_out, t_out, r_out ] = BiRed( A, t, r )

  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );

  [ tT, ...
    tB ] = FLA_Part_2x1( t, ...
                         0, 'FLA_TOP' );

  [ rT, ...
    rB ] = FLA_Part_2x1( r, ...
                         0, 'FLA_TOP' );
                     
  while ( size( ATL, 1 ) < size( A, 1 ) )

    [ A00,  a01,     A02,  ...
      a10t, alpha11, a12t, ...
      A20,  a21,     A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                                    ABL, ABR, ...
                                                    1, 1, 'FLA_BR' );

    [ t0, ...
      tau1, ...
      t2 ] = FLA_Repart_2x1_to_3x1( tT, ...
                                    tB, ...
                                    1, 'FLA_BOTTOM' );

    [ r0, ...
      rho1, ...
      r2 ] = FLA_Repart_2x1_to_3x1( rT, ...
                                    rB, ...
                                    1, 'FLA_BOTTOM' );
                                
    %------------------------------------------------------------%
    % Step 1: Construct the column vector x_col by combining the current 
    % diagonal element alpha11 and the subdiagonal elements a21.
    x_col = [alpha11; a21];
    
    % Step 2: Compute the Householder vector u_col and the scalar tau1
    % using the Housev1 function. 
    % - u_col is the Householder vector that will be used to zero out the 
    %   elements of x_col (except for the first element).
    % - tau1 is a scalar used in the transformation.
    [u_col, tau1] = Housev1(x_col);
    
    % Store tau1 in the t vector
    t(size(tT, 1) + 1) = tau1;
    
    % Step 3: Update a21 and A22 using the Householder vector
    v = u_col(2:end);
    alpha11 = u_col(1);  % Update the diagonal element
    a21 = v;
    
    if size(A22, 1) > 1 && size(A22, 2) > 1
        % Form the full Householder vector u_col
        u_col_full = [1; v];
        
        % Apply the Householder transformation to [a12t; A22]
        disp(size(a12t'))
        disp(size(A22))
        block_matrix = [a12t; A22'];
        
        % Householder matrix H = I - (1/tau1) * u_col_full * u_col_full'
        H_col = eye(size(block_matrix, 1)) - (1/tau1) * (u_col_full * u_col_full');
        
        % Apply the transformation to the block matrix
        transformed_block = H_col * block_matrix;
        
        % Extract the updated row vector a12t and submatrix A22
        a12t = transformed_block(1, :);  % Updated row vector
        A22 = transformed_block(2:end, :);  % Updated submatrix
    end


    
    % Zero out all of a21
    a21(:) = 0;
    
    
    % Reassemble the matrix before proceeding to the row transformation
    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
                                             a10t, alpha11, a12t, ...
                                             A20,  a21,     A22, ...
                                             'FLA_TL' );

    % Everything Until Here Looks good
    
    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
                                             a10t, alpha11, a12t, ...
                                             A20,  a21,     A22, ...
                                             'FLA_TL' );

    [ tT, ...
      tB ] = FLA_Cont_with_3x1_to_2x1( t0, ...
                                       tau1, ...
                                       t2, ...
                                       'FLA_TOP' );

    [ rT, ...
      rB ] = FLA_Cont_with_3x1_to_2x1( r0, ...
                                       rho1, ...
                                       r2, ...
                                       'FLA_TOP' );
                                   
  end

  A_out = [ ATL, ATR
            ABL, ABR ];

  t_out = [ tT
            tB ];
        
  r_out = [ rT
            rB ];

return