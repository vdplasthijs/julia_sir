### A Pluto.jl notebook ###
# v0.7.5

using Markdown
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.peek, el) ? Base.peek(el) : missing
        el
    end
end
# ╔═╡ 0b10b4d6-8989-11ea-0b19-e35480058578
begin 
	n_timepoints = 500
	time_first_inf = 30
	n_first_inf = 20
	n_pop = 10000
	time_recovery = 1
	""
end

# ╔═╡ dc771cce-8993-11ea-0810-497feac7d2f6
begin
	function forward_outbreak(; r_inf=0.05, r_dec=0.01, r_rec=0.01, r_vac=0.1)
		susc = zeros(n_timepoints)
		inf = zeros(n_timepoints)
		rec = zeros(n_timepoints)
		dec= zeros(n_timepoints)
		susc[1:time_first_inf] = ones(time_first_inf) * n_pop * (1 - r_vac)
		rec[1:time_first_inf] = ones(time_first_inf) * n_pop * r_vac
		# model 
		for t in (time_recovery + 1):n_timepoints
			susc[t] = susc[t - 1] - r_inf * inf[t - 1] * susc[t - 1] / n_pop
			if t == time_first_inf
				inf[t] = n_first_inf
			else
				inf[t] = inf[t - 1] + r_inf * inf[t - 1] * susc[t - 1] / n_pop - r_rec * inf[t - time_recovery] - r_dec * inf[t - time_recovery]
			end
			rec[t] = rec[t - 1] + r_rec * inf[t - time_recovery]
			dec[t] = dec[t - 1] + r_dec * inf[t - time_recovery]
		end
		return (susc, inf, rec, dec)
	end
	md"## SIR Model
	
	Very interesting stuff!! Check colourful plots below"
end

# ╔═╡ a23510e8-89f5-11ea-0c25-a75768612b3b
begin
	inf_slider = @bind rate_inf html"<input type='range' min='0.0' max='1' step='0.01' value='0.15'>"
	rec_slider = @bind rate_rec html"<input type='range' min='0.0' max='1' step='0.01' value='0.04'>"
	dec_slider = @bind rate_dec html"<input type='range' min='0.0' max='1' step='0.01' value='0.02'>"
	vac_slider = @bind rate_vac html"<input type='range' min='0.0' max='1' step='0.01' value='0.0'>"
	
	md"""**Please set model parameters:**
	
	Infection rate: $(inf_slider)
	
	Recovery rate: $rec_slider
	
	Decease rate: $(dec_slider)

	Vaccination rate: $(vac_slider)
	"""
end

# ╔═╡ c1095b84-8a20-11ea-1fe4-f5d22fc9db85

begin
	
	survival = rate_rec / (rate_rec + rate_dec)
	reproduction = rate_inf / (rate_rec + rate_dec)
	
	md"""Currently $R_0$ = $(round(reproduction, digits=3)), and survival fraction is $(round(survival, digits=3)). Change parametesr below to alter"""
end

# ╔═╡ a731f220-89f8-11ea-311e-01591861953a
md"""Infection rate: $rate_inf, Recovery rate: $rate_rec, 
Decease rate: $rate_dec, Vaccination rate: $rate_vac"""

# ╔═╡ 65d2d834-898b-11ea-291c-f3ccc5317e69
begin
	ex_susc, ex_inf, ex_rec, ex_dec = forward_outbreak(; r_inf=rate_inf, r_dec=rate_dec, r_rec=rate_rec, r_vac=rate_vac)
	floor_plot = ones(n_timepoints) * 0
	ceiling_plot = ones(n_timepoints) * n_pop
	plot(ex_inf + ex_susc, ribbon=(ex_susc, floor_plot), label="Susceptible")
	plot!(ex_inf, ribbon=(ex_inf, floor_plot), xlabel="Time", 
		ylabel="Population", label="Infected", title="Epidemic evolution")
	plot!(ex_inf + ex_susc + ex_rec, ribbon=(ex_rec, floor_plot), label="Recovered")
	plot!(ceiling_plot, ribbon=(ex_dec, floor_plot), label="Deceased")
end

# ╔═╡ 75e9bc76-898e-11ea-38ed-6591b6b0ff69
begin
	vac_array = 0:0.01:1
	dec_array = zeros(length(vac_array))
	for i_v in 1:length(vac_array)
		s, i, r, d = forward_outbreak(; r_inf=rate_inf, r_rec=rate_rec, r_dec=rate_dec, r_vac=vac_array[i_v])  # run simulation
		dec_array[i_v] = d[end] / n_pop # final number of deaths
	end
	
	plot(vac_array, dec_array, xlabel="Fraction vaccinated", 
		ylabel="Fraction deceased", lw=3, title="Number of deceased depends nonlinearly on vaccination", label="Deceased (fraction)")
end

# ╔═╡ f1e677be-8a20-11ea-0c03-3bd72078ae95
begin
	rate_dec_array = 0:0.005:1
	dec_array_1 = zeros(length(rate_dec_array))
	for i_r in 1:length(rate_dec_array)
		s1, i1, r1, d1 = forward_outbreak(; r_inf=rate_inf, r_rec=rate_rec, r_dec=rate_dec_array[i_r], r_vac=rate_vac)  # run simulation
		dec_array_1[i_r] = d1[end] / n_pop # final number of deaths
	end
	
	plot(rate_dec_array, dec_array_1, xlabel="Mortality rate", 
		ylabel="Fraction deceased", lw=3, title="Number of deceased depends non-monotonically on \nmortality rate", label="Deceased (fraction)")
end

# ╔═╡ 35568376-898e-11ea-2293-7de50f09cea5
begin
	using Plots
	md""" # Bonjour"""
end

# ╔═╡ Cell order:
# ╟─35568376-898e-11ea-2293-7de50f09cea5
# ╟─0b10b4d6-8989-11ea-0b19-e35480058578
# ╟─c1095b84-8a20-11ea-1fe4-f5d22fc9db85
# ╟─dc771cce-8993-11ea-0810-497feac7d2f6
# ╟─a23510e8-89f5-11ea-0c25-a75768612b3b
# ╟─a731f220-89f8-11ea-311e-01591861953a
# ╟─65d2d834-898b-11ea-291c-f3ccc5317e69
# ╟─75e9bc76-898e-11ea-38ed-6591b6b0ff69
# ╟─f1e677be-8a20-11ea-0c03-3bd72078ae95
