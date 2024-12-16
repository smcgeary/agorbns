import statsmodels.api as sm
import statsmodels.formula.api as smf


data = sm.datasets.get_rdataset("dietox", "geepack").data
print(data)
md = smf.mixedlm("Weight ~ Time", data, groups=data["Pig"])

mdf = md.fit()

print(mdf.summary())