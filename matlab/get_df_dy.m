function df_dy = get_df_dy(t,x,y,ps)

[~,~,df_dy] = differential_eqs(t,x,y,ps);

