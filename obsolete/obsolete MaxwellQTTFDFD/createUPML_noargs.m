%createUPML_noargs
if doUPML
    
    
    sp_x=sparse(full(syz_x));
    sp_y=sparse(full(sxz_y));
    sp_z=sparse(full(sxy_z));

    eps_Ex_f=sparse(eps_3_L.*sp_x);
    eps_Ey_f=sparse(eps_3_L.*sp_y);
    eps_Ez_f=sparse(eps_3_L.*sp_z);

    eps_Ex_m_f=sparse(diag(eps_Ex_f));
    eps_Ey_m_f=sparse(diag(eps_Ey_f));
    eps_Ez_m_f=sparse(diag(eps_Ez_f));

    eps_All_f=[      eps_Ex_m_f,    nil_f,        nil_f;
                     nil_f,         eps_Ey_m_f,   nil_f;
                     nil_f,         nil_f,        eps_Ez_m_f     ];

    invspyz_x=sparse(full(invsyz_x));
    invspxz_y=sparse(full(invsxz_y));
    invspxy_z=sparse(full(invsxy_z));

    inv_eps_3_L=mask_of_coords_f*eps^(-1) + (1-mask_of_coords_f);

    inv_eps_Ex_f=sparse(inv_eps_3_L.*invspyz_x);
    inv_eps_Ey_f=sparse(inv_eps_3_L.*invspxz_y);
    inv_eps_Ez_f=sparse(inv_eps_3_L.*invspxy_z);

    inv_eps_Ex_m_f=diag(inv_eps_Ex_f);
    inv_eps_Ey_m_f=diag(inv_eps_Ey_f);
    inv_eps_Ez_m_f=diag(inv_eps_Ez_f);

    inv_eps_All_f=[inv_eps_Ex_m_f, nil_f,          nil_f;
                nil_f,         inv_eps_Ey_m_f,  nil_f;
                nil_f,         nil_f,          inv_eps_Ez_m_f     ];


    inv_mu_Ex_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));
    inv_mu_Ey_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));
    inv_mu_Ez_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));

    inv_mu_Ex_m_f=diag(inv_mu_Ex_f.*invspyz_x);
    inv_mu_Ey_m_f=diag(inv_mu_Ey_f.*invspxz_y);
    inv_mu_Ez_m_f=diag(inv_mu_Ez_f.*invspxy_z);


    inv_mu_All_f=[inv_mu_Ex_m_f, nil_f,        nil_f;
                nil_f,         inv_mu_Ey_m_f,  nil_f;
                nil_f,         nil_f,          inv_mu_Ez_m_f     ];
else
    eps_Ex_f=sparse(eps_3_L);
    eps_Ey_f=sparse(eps_3_L);
    eps_Ez_f=sparse(eps_3_L);

    eps_Ex_m_f=sparse(diag(eps_Ex_f));
    eps_Ey_m_f=sparse(diag(eps_Ey_f));
    eps_Ez_m_f=sparse(diag(eps_Ez_f));

    eps_All_f=[      eps_Ex_m_f,    nil_f,        nil_f;
                     nil_f,         eps_Ey_m_f,   nil_f;
                     nil_f,         nil_f,        eps_Ez_m_f     ];

    inv_mu_Ex_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));
    inv_mu_Ey_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));
    inv_mu_Ez_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));

    inv_mu_Ex_m_f=diag(inv_mu_Ex_f);
    inv_mu_Ey_m_f=diag(inv_mu_Ey_f);
    inv_mu_Ez_m_f=diag(inv_mu_Ez_f);

    inv_mu_All_f=[inv_mu_Ex_m_f, nil_f,        nil_f;
                nil_f,         inv_mu_Ey_m_f,  nil_f;
                nil_f,         nil_f,          inv_mu_Ez_m_f     ];             

end