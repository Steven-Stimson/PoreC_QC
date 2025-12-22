#!/usr/bin/env perl
use strict;
use warnings;

# 默认参数
my $add_commas = 0;
my $decimal_places = -1;

# 解析命令行参数
while (@ARGV && $ARGV[0] =~ /^-/) {
    my $arg = shift @ARGV;
    if ($arg eq '-k') {
        $add_commas = 1;
    } elsif ($arg eq '-f' && @ARGV) {
        $decimal_places = shift @ARGV;
        die "Invalid number of decimal places: $decimal_places\n" unless $decimal_places =~ /^\d+$/;
    } else {
        die "Unknown option: $arg\n";
    }
}

# 检查标准输入是否有内容
my $input = do { local $/; <STDIN> };  # 读取全部标准输入内容
if (defined $input && $input ne '') {
    # 标准输入有内容，直接处理
    process_input($input);
} elsif (@ARGV) {
    # 标准输入为空，尝试打开第一个参数作为文件
    my $file = shift @ARGV;
    open my $fh, '<', $file or die "Cannot open $file: $!";
    my $file_content = do { local $/; <$fh> };
    close $fh;
    process_input($file_content);
} else {
    die "No input provided via stdin or file argument.\n";
}

# 处理输入内容的子函数
sub process_input {
    my ($content) = @_;
    foreach my $line (split /\n/, $content) {
        chomp $line;
        my @fields = split /\t/, $line;  # 按照 tab 分割

        # 遍历每个字段进行处理
        foreach my $field (@fields) {
            if ($field =~ /^\d+$/) {  # 如果字段是纯数字
                $field = digit_formatter($field) if $add_commas;
            } elsif ($field =~ /^-?\d+\.\d+$/) {  # 如果字段是浮点数
                if ($decimal_places >= 0) {
                    $field = sprintf("%.${decimal_places}f", $field);  # 四舍五入到指定小数位数
                }
                $field = digit_formatter($field) if $add_commas;
            }
        }

        # 使用 join 将处理后的字段用 Tab 连接，并输出
        print join("\t", @fields)."\n";
    }
}

# 数字格式化函数
sub digit_formatter {
    my ($num) = @_;
    if ($num =~ /\./) {  # 如果是浮点数
        my ($integer, $decimal) = split /\./, $num;
        $integer = reverse $integer;
        $integer =~ s/(...)/$1,/g;
        $integer = reverse $integer;
        $integer =~ s/^,//;
        return "$integer.$decimal";
    } else {  # 如果是整数
        $num = reverse $num;
        $num =~ s/(...)/$1,/g;
        $num = reverse $num;
        $num =~ s/^,//;
        return $num;
    }
}
