function Code(el)
    return pandoc.RawInline('latex', '\\inlinecode{' .. el.text .. '}')
  end
