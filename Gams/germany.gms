* Define a set for our time step in our model.
* Our model has an hourly resolution, and therefore 8760 time steps. 
set h/
    h1*h8760
/;

* Define a set for our technologies, t
set t/
    wind,
    pv,
    gas,
    gas_ccs
/;

* Define a set for our storage technologies, s
set s/
    battery
/;

parameter demand(h)/
* hourly demand for new capacity in Germany in MWh
* calculated by scaling the demand with the fraction of the electricity that needs to be phased out (coal and nuclear fraction of electricity mix)
$include Germany_altdemand_2019.tsv
/;

parameter capacity_factors(h,t)/
* capacity factors for pv, wind, gas and gas_ccs
$include germany_capacity_factors_2019.tsv
/;

parameter annuitised_costs(t)/
* unit is EUR/MW
* discount rate used is 0.05
* Source for gas and gas_ccs current cost: https://www.globalccsinstitute.com/archive/hub/publications/17011/costs-co2-capture-transport-and-storage.pdf
* Source for gas and gas_ccs cost forecast to 2050: https://www.sciencedirect.com/science/article/pii/S2211467X18300634?via%3Dihub#tbl14
* Source for wind and PV cost: https://www.iea.org/reports/world-energy-model/techno-economic-inputs

* Annuitised costs
* year         2020,  2025,  2030,  2035,  2040,  2045,  2050

* Wind:       85868, 84628, 83387, 82147, 80907, 79666, 78426
* PV:         52448, 48181, 43915, 39648, 35382, 31115, 26848
* CCGT:       44995, 43956, 42918, 41880, 40842, 39803, 38765
* CCGT + CCS: 86498, 83389, 80280, 77172, 74063, 70955, 67846

wind 78426
pv 26848
gas 38765
gas_ccs 67846
/;

parameter variable_costs(t)/
* EUR/MWh
* fuel cost for gas: https://github.com/PyPSA/pypsa-eur/blob/master/data/costs.csv
* Efficiencies for gas and gas_ccs: https://www.globalccsinstitute.com/archive/hub/publications/201688/global-ccs-cost-updatev4.pdf
* var_cost = (fuel_cost / efficiency) + VOM      ( (EUR/MWh_th) / (MWh_el/MWh_th) = EUR/MWh_el )
wind 2.3
pv 0.01
gas 45.9
gas_ccs 51.3
/;

parameter availability(t)/
* Source: https://energymag.net/technologies/availability-factor/
wind 0.95
pv 0.98
gas 0.9
gas_ccs 0.9
/;

parameter emission_factor(t)/
* tCO/MWh
* from: https://www.globalccsinstitute.com/archive/hub/publications/201688/global-ccs-cost-updatev4.pdf
wind 0
pv 0
gas 0.356
gas_ccs 0.04
/;

parameter storage_cost(s)/
* Annuitised costs for the storage unit [EUR/MWh]
* Source for investment costs: https://www.nrel.gov/docs/fy21osti/79236.pdf
* Costs 2020, 2025, 2030, 2035, 2040, 2045, 2050
* 29250, 20517, 16787, 15769, 14752, 13650, 12632
battery 12632
/;

scalar
emission_cost
* EUR/tCO2, based on ETS price of CO2
* Prices with current policies: 48.4, 52.8, 57.2, 61.6, 66, 73, 79.2
* Prcies in net-zero scenario: 48.4, 81.84, 114.4, 146.96, 180.4,  200.64, 220
/79.2/
;

positive variables
    var_installed_cap(t)
    var_g(h,t)
    var_emissions(t)
    var_invest_cost(t)
    var_var_cost(t)
    var_emission_cost(t)
* Storage variables
    var_installed_storage(s)
    var_storage_level(h,s)
    var_storage_cost(s)
    var_charge(h,s)
    var_discharge(h,s)
;

free variables
* default so we could just write "variables"
* the thing we maximise/minimise need to be a free variable (although it is a positive variable). 
    var_system_cost
; 

equations
    eq_objective
    eq_demand_balance
    eq_gen_convert
    eq_investment_cost
    eq_variable_costs
    eq_emission_costs
    eq_emissions
* Storage equations
    eq_storage_level
    eq_storage_cost
    eq_storage_constraint
;

*** OBJECTIVE FUNCTION ***
eq_objective.. var_system_cost =E=
* Investment costs
    sum((t),annuitised_costs(t) * var_installed_cap(t)) +
* Variable costs
    sum((h,t),variable_costs(t) * var_g(h,t)) +
* Storage costs
    sum((s),storage_cost(s) * var_installed_storage(s)) +
* emission costs
    sum((t), var_emission_cost(t))
;

* INVESTMENT COSTS
eq_investment_cost(t).. var_invest_cost(t) =E= var_installed_cap(t)*annuitised_costs(t);

* VARIABLE COSTS
eq_variable_costs(t).. var_var_cost(t) =E= sum(h,variable_costs(t)*var_g(h,t));

* RELATIONSHIP BETWEEN CAPACITY FACTORS AND GENERATION
eq_gen_convert(h,t).. var_g(h,t) =L= var_installed_cap(t) * capacity_factors(h,t) * availability(t);

* STORAGE COSTS
eq_storage_cost(s).. var_storage_cost(s) =E= var_installed_storage(s)*storage_cost(s);

* DEMAND BALANCE EQUATION
eq_demand_balance(h).. sum((t),var_g(h,t)) + sum((s),var_discharge(h,s)) - sum((s),var_charge(h,s)) =E= demand(h);

* CO2 EMITTED
eq_emissions(t).. var_emissions(t) =E= sum((h),var_g(h,t)*emission_factor(t));
* MWh * tCO2e/MWh = tCO2e

* EMISSION COSTS
eq_emission_costs(t).. var_emission_cost(t) =E= emission_cost * sum((h),var_g(h,t)*emission_factor(t));

* CURRENT STORAGE LEVEL
eq_storage_level(h,s).. var_storage_level(h,s) =E=
    var_storage_level(h-1,s) +
* balance between charge and discharge
    var_charge(h,s) - var_discharge(h,s);

* STORAGE CONSTRAINT
eq_storage_constraint(h,s).. var_storage_level(h,s) =L= var_installed_storage(s);

model optimal_generation /
                            all/;
                            
solve optimal_generation minimizing var_system_cost using lp;

execute_unload 'germany.gdx' var_g, var_installed_cap, var_storage_level, var_installed_storage, var_charge, var_discharge, var_emissions, var_system_cost;
execute "gdx2sqlite -i germany.gdx -o model_results2050.db"