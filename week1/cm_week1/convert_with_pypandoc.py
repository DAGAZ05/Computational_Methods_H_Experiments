import pypandoc

# 最简单的用法：文件到文件
def simple_convert():
    # 将 markdown 转换为 word
    output = pypandoc.convert_file('test.md', 'docx', outputfile='test.docx')
    print("转换完成！")


# 转换字符串到文件
def string_to_file():
    md_text = """
    # 标题

    这是一个 **测试**。

    $\\alpha + \\beta = \\gamma$
    """

    pypandoc.convert_text(md_text, 'docx', format='md', outputfile='output.docx')


# 获取转换后的内容（不进文件）
def get_converted_content():
    md_text = "# Hello *World*"

    # 转换为 HTML
    html = pypandoc.convert_text(md_text, 'html', format='md')
    print(html)  # <h1>Hello <em>World</em></h1>

    # 转换为 Word 的二进制内容
    docx_bytes = pypandoc.convert_text(md_text, 'docx', format='md')
    # 可以保存或进一步处理
    with open('output.docx', 'wb') as f:
        f.write(docx_bytes)


# 使用模板
def convert_with_template():
    pypandoc.convert_file(
        'input.md',
        'docx',
        outputfile='output.docx',
        extra_args=['--reference-doc=template.docx']
    )


# 高级选项
def advanced_convert():
    # 添加多个参数
    output = pypandoc.convert_file(
        'report.md',
        'docx',
        outputfile='report.docx',
        extra_args=[
            # '--toc',  # 生成目录
            '--number-sections',  # 章节编号
            '--mathml',  # 使用 MathML 处理公式
            # '-V', 'geometry:margin=1in'  # 设置页边距
        ]
    )

if __name__ == "__main__":
    simple_convert()
    # advanced_convert()
