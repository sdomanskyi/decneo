from scRegulation.commonFunctions import *

if __name__ == '__main__':

    wdir = 'd:/Projects/A_Endothelial/VS/Endothelial/results/c/'

    dfD = pd.read_excel(wdir + 'summary all PanglaoDB.xlsx', sheet_name='DCS', inde_col=None, header=0)
    dfA = pd.read_excel(wdir + 'summary all PanglaoDB.xlsx', sheet_name='Alona', inde_col=None, header=0)

    def getM(df, s):
        df = pd.concat([df[df.columns[0:2]].set_index('3m').dropna(), 
                       df[df.columns[2:4]].set_index('4m').dropna(), 
                       df[df.columns[4:6]].set_index('3h').dropna(), 
                       df[df.columns[6:]].set_index('4h').dropna()], 
                       axis=1, sort=False, keys=['3m' + s, '4m' + s, '3h' + s, '4h' + s]).fillna(0.).droplevel(1, axis=1)

        return df

    dfD = getM(dfD, '_DCS')
    dfA = getM(dfA, '_Alona')
    df = pd.concat([dfA, dfD], axis=1, sort=False).fillna(0.)
    print(df)

    df.to_excel(wdir + 'summary_aligned.xlsx')
    
    exit()




    print(dfA)

