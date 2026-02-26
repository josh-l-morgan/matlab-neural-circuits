MyCZEMAPIClass = actxserver('VBComObjectWrapperForZeissAPI.KHZeissSEMWrapperComClass')

MyCZEMAPIClass.InitialiseRemoting


MyCZEMAPIClass.GetMag
MyCZEMAPIClass.Get_ReturnTypeSingle('AP_MAG')
MyCZEMAPIClass.Get_ReturnTypeString('DP_COLUMN_TYPE')

MyCZEMAPIClass.Set_PassedTypeSingle('AP_MAG',100)
