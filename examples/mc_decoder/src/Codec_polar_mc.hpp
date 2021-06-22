/*!
 * \file
 * \brief Class factory::Codec_polar_mc.
 */
// 宏定义，控制代码是否编译（预处理器指令中符号常量命名通常约定头文件名用大写形式，圆点用下划线代替
// FACTORY::CODEC_POLAR_MC.HPP
#ifndef FACTORY_CODEC_POLAR_MC_HPP
#define FACTORY_CODEC_POLAR_MC_HPP

#include <vector>
#include <string>
#include <map>
#include <cli.hpp>

#include <aff3ct.hpp>


//#include "Tools/Factory/Header.hpp"
//#include "Module/CRC/CRC.hpp"
//#include "Tools/Codec/Polar/Codec_polar.hpp"
//#include "Factory/Tools/Code/Polar/Frozenbits_generator.hpp"

// 包含了系统的基类Codec_SISO
//#include "Factory/Tools/Codec/Codec_SISO.hpp"

// 名称空间知识点
namespace aff3ct
{
namespace factory
{
// 声明与定义的区别,此处是声明
extern const std::string Codec_polar_name;
extern const std::string Codec_polar_prefix;

// Codec_polar 由 Codec_SISO 继承而来。涉及到类的继承
// Codec_SISO是基类，public表示共有派生，基类的公有成员将成为派生类的公有成员，基类的私有部分只能通过基类的公有和保护方法访问
class Codec_polar_mc : public Codec_SISO
{
public:
	// ----------------------------------------------------------------------------------------------------- PARAMETERS
	// depending parameters
	// 智能指针，自动销毁
	tools::auto_cloned_unique_ptr<Frozenbits_generator> fbg;

	// -------------------------------------------------------------------------------------------------------- METHODS
	// explicit：指定构造函数或转换函数 (C++11起)为显式, 即它不能用于隐式转换和复制初始化.
	explicit Codec_polar_mc(const std::string &p = Codec_polar_prefix);
	// Virtual关键字的函数就是虚拟函数，于是在Base的派生类Derived中就可以通过重写虚拟函数来实现对基类虚拟函数的覆盖。
	virtual ~Codec_polar_mc() = default;
	Codec_polar_mc* clone() const;
	void enable_puncturer();

	virtual std::vector<std::string> get_names      () const;
	virtual std::vector<std::string> get_short_names() const;
	virtual std::vector<std::string> get_prefixes   () const;

	// parameters construction
	void get_description(cli::Argument_map_info &args) const;
	void store          (const cli::Argument_map_value &vals);
	void get_headers    (std::map<std::string,tools::header_list>& headers, const bool full = true) const;

	// template 函数模板
	// builder
	template <typename B = int, typename Q = float>
	tools::Codec_polar<B,Q>* build(const module::CRC<B>* crc = nullptr) const;
};
}
}

#endif /* FACTORY_CODEC_POLAR_MC_HPP */
