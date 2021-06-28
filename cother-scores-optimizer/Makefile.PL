use ExtUtils::MakeMaker;
WriteMakefile(
    "NAME" => "ddms",#package name
    "CC" => "g++",
    "LD" => "g++",
    "OPTIMIZE" => "-O3",
    "OBJECT" => "ddms.o ddms_wrap.o"
);
