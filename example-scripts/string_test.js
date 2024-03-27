
print(format("format(str,...): str='{0}'","str"));
print(format("format(str,...): str='{0}' int(123)='{1}' float(123.456)='{2}'","str", 123, 123.456));
print("str.format(...): str='{0}' int(123)='{1}' float(123.456)='{2}'".format("str", 123, 123.456));
var sfmt = "%03i";
print("fmti({0},12) = {1}".format(sfmt,fmti(sfmt,12)));
print("fmti('abi_A%03i_B%03i,12) = '{0}'".format(fmti("abi_A%03i_B%03i",12)));
print("fmti(fmti('abi_A%03i_B%03i,12),34) = '{0}'".format(fmti(fmti("abi_A%03i_B%03i",12),34)));
// ^-- will produce garbage if called with more than one argument..
